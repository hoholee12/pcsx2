
/*------------------------------------------------------------------------------
                            [GLOBALS|FUNCTIONS]
------------------------------------------------------------------------------*/

float2 texcoord = GetCoordinates();
float2 pixelSize = GetInvResolution();
float2 screenSize = GetResolution();
float2 texSize = textureSize(samp0, 0).xy;

#define float3x3 mat3
#define mul(x, y) y * x
#define FIX(c) max(abs(c), 1e-5)
#define saturate(x) clamp(x, 0.0, 1.0)
#define SampleLocationLod(location, lod) textureLod(samp0, float3(location, layer), lod)

float Epsilon = 1e-10;
const float3 lumCoeff = float3(0.2126729, 0.7151522, 0.0721750);

//Average relative luminance
float AvgLuminance(float3 color)
{
    return sqrt(
    (color.x * color.x * lumCoeff.x) +
    (color.y * color.y * lumCoeff.y) +
    (color.z * color.z * lumCoeff.z));
}

float smootherstep(float a, float b, float x)
{
    x = saturate((x - a) / (b - a));
    return x*x*x*(x*(x * 6.0 - 15.0) + 10.0);
}


/*------------------------------------------------------------------------------
                        [BICUBIC SCALER CODE SECTION]
------------------------------------------------------------------------------*/

float4 BicubicScaler(float2 uv, float2 txSize)
{
    float2 inputSize = float2(1.0/txSize.x, 1.0/txSize.y);

    float2 coord_hg = uv * txSize - 0.5;
    float2 index = floor(coord_hg);
    float2 f = coord_hg - index;

    mat4 M = mat4( -1.0, 3.0,-3.0, 1.0, 3.0,-6.0, 3.0, 0.0,
                   -3.0, 0.0, 3.0, 0.0, 1.0, 4.0, 1.0, 0.0 );
    M /= 6.0;

    float4 wx = M * float4(f.x*f.x*f.x, f.x*f.x, f.x, 1.0);
    float4 wy = M * float4(f.y*f.y*f.y, f.y*f.y, f.y, 1.0);
    float2 w0 = float2(wx.x, wy.x);
    float2 w1 = float2(wx.y, wy.y);
    float2 w2 = float2(wx.z, wy.z);
    float2 w3 = float2(wx.w, wy.w);

    float2 g0 = w0 + w1;
    float2 g1 = w2 + w3;
    float2 h0 = w1 / g0 - 1.0;
    float2 h1 = w3 / g1 + 1.0;

    float2 coord00 = index + h0;
    float2 coord10 = index + float2(h1.x, h0.y);
    float2 coord01 = index + float2(h0.x, h1.y);
    float2 coord11 = index + h1;

    coord00 = (coord00 + 0.5) * inputSize;
    coord10 = (coord10 + 0.5) * inputSize;
    coord01 = (coord01 + 0.5) * inputSize;
    coord11 = (coord11 + 0.5) * inputSize;

    float4 tex00 = SampleLocationLod(coord00, 0.0);
    float4 tex10 = SampleLocationLod(coord10, 0.0);
    float4 tex01 = SampleLocationLod(coord01, 0.0);
    float4 tex11 = SampleLocationLod(coord11, 0.0);

    tex00 = lerp(tex01, tex00, float4(g0.y, g0.y, g0.y, g0.y));
    tex10 = lerp(tex11, tex10, float4(g0.y, g0.y, g0.y, g0.y));

    float4 res = lerp(tex10, tex00, float4(g0.x, g0.x, g0.x, g0.x));

    return res;
}

float4 BicubicScalerPass(float4 color)
{
    color = BicubicScaler(texcoord, screenSize);
    return color;
}

/*------------------------------------------------------------------------------
                         [LANCZOS SCALER CODE SECTION]
------------------------------------------------------------------------------*/

float3 PixelPos(float xpos, float ypos)
{
    return SampleLocation(float2(xpos, ypos)).rgb;
}

float4 WeightQuad(float x)
{
    const float PI = 3.141592653;
    float4 weight = FIX(PI * float4(1.0 + x, x, 1.0 - x, 2.0 - x));
    float4 ret = sin(weight) * sin(weight / 2.0) / (weight * weight);

    return ret / dot(ret, float4(1.0, 1.0, 1.0, 1.0));
}

float3 LineRun(float ypos, float4 xpos, float4 linetaps)
{
    return mat4x3(
    PixelPos(xpos.x, ypos),
    PixelPos(xpos.y, ypos),
    PixelPos(xpos.z, ypos),
    PixelPos(xpos.w, ypos)) * linetaps;
}

float4 LanczosScaler(float2 inputSize)
{
    float2 stepxy = float2(1.0/inputSize.x, 1.0/inputSize.y);
    float2 pos = texcoord + stepxy * 0.5;
    float2 f = frac(pos / stepxy);

    float2 xystart = (-1.5 -f) * stepxy + pos;
    float4 xpos = float4(
    xystart.x,
    xystart.x + stepxy.x,
    xystart.x + stepxy.x * 2.0,
    xystart.x + stepxy.x * 3.0);

    float4 linetaps = WeightQuad(f.x);
    float4 columntaps = WeightQuad(f.y);

    return float4(mat4x3(
    LineRun(xystart.y, xpos, linetaps),
    LineRun(xystart.y + stepxy.y, xpos, linetaps),
    LineRun(xystart.y + stepxy.y * 2.0, xpos, linetaps),
    LineRun(xystart.y + stepxy.y * 3.0, xpos, linetaps)) * columntaps, 1.0);
}

float4 LanczosScalerPass(float4 color)
{
    color = LanczosScaler(screenSize);
    return color;
}

/*------------------------------------------------------------------------------
                       [BILINEAR FILTERING CODE SECTION]
------------------------------------------------------------------------------*/

float4 SampleBilinear(float2 texcoord)
{
    float texelSizeX = 1.0 / texSize.x;
    float texelSizeY = 1.0 / texSize.y;

    int nX = int(texcoord.x * texSize.x);
    int nY = int(texcoord.y * texSize.y);

    float2 uvCoord = float2(float(nX) / texSize.x, float(nY) / texSize.y);

    // Take nearest two data in current row.
    float4 SampleA = SampleLocation(uvCoord);
    float4 SampleB = SampleLocation(uvCoord + float2(texelSizeX, 0.0));

    // Take nearest two data in bottom row.
    float4 SampleC = SampleLocation(uvCoord + float2(0.0, texelSizeY));
    float4 SampleD = SampleLocation(uvCoord + float2(texelSizeX , texelSizeY));

    float LX = fract(texcoord.x * texSize.x); //interpolation factor for X direction.
    float LY = fract(texcoord.y * texSize.y); //interpolation factor for Y direction.

    // Interpolate in X direction.
    float4 InterpolateA = mix(SampleA, SampleB, LX); //Top row in X direction.
    float4 InterpolateB = mix(SampleC, SampleD, LX); //Bottom row in X direction.

    return mix(InterpolateA, InterpolateB, LY); //Interpolate in Y direction.
}

float4 BilinearPass(float4 color)
{
    color = SampleBilinear(texcoord);
    return color;
}


/*------------------------------------------------------------------------------
                       [PIXEL VIBRANCE CODE SECTION]
------------------------------------------------------------------------------*/

float A_VIBRANCE = 0.15;
float B_R_VIBRANCE = 1.00;
float C_G_VIBRANCE = 1.00;
float D_B_VIBRANCE = 1.00;

float4 VibrancePass(float4 color)
{
    float vib = A_VIBRANCE;
    float luma = AvgLuminance(color.rgb);

    float colorMax = max(color.r, max(color.g, color.b));
    float colorMin = min(color.r, min(color.g, color.b));

    float colorSaturation = colorMax - colorMin;
    float3 colorCoeff = float3(B_R_VIBRANCE*
    vib, C_G_VIBRANCE * vib, D_B_VIBRANCE * vib);

    color.rgb = lerp(float3(luma), color.rgb, (1.0 + (colorCoeff * (1.0 - (sign(colorCoeff) * colorSaturation)))));
    color.a = AvgLuminance(color.rgb);

    return saturate(color); //Debug: return colorSaturation.xxxx;
}


/*------------------------------------------------------------------------------
                       [TEXTURE SHARPEN CODE SECTION]
------------------------------------------------------------------------------*/


float A_SHARPEN_STRENGTH = 0.75;
float B_SHARPEN_CLAMP = 0.015;
float C_SHARPEN_BIAS = 1.20;


float Cubic(float coeff)
{
    float4 n = float4(1.0, 2.0, 3.0, 4.0) - coeff;
    float4 s = n * n * n;

    float x = s.x;
    float y = s.y - 4.0 * s.x;
    float z = s.z - 4.0 * s.y + 6.0 * s.x;
    float w = 6.0 - x - y - z;

    return (x + y + z + w) / 4.0;
}

float4 SampleBicubic(float2 texcoord)
{
    float texelSizeX = (1.0 / texSize.x) * C_SHARPEN_BIAS;
    float texelSizeY = (1.0 / texSize.y) * C_SHARPEN_BIAS;

    float4 nSum = float4(0.0, 0.0, 0.0, 0.0);
    float4 nDenom = float4(0.0, 0.0, 0.0, 0.0);

    float a = frac(texcoord.x * texSize.x);
    float b = frac(texcoord.y * texSize.y);

    int nX = int(texcoord.x * texSize.x);
    int nY = int(texcoord.y * texSize.y);

    float2 uvCoord = float2(float(nX) / texSize.x, float(nY) / texSize.y);

    for (int m = -1; m <= 2; m++) {
    for (int n = -1; n <= 2; n++) {
    
    float4 Samples = SampleLocation(uvCoord +
    float2(texelSizeX * float(m), texelSizeY * float(n)));

    float vc1 = Cubic(float(m) - a);
    float4 vecCoeff1 = float4(vc1, vc1, vc1, vc1);

    float vc2 = Cubic(-(float(n) - b));
    float4 vecCoeff2 = float4(vc2, vc2, vc2, vc2);

    nSum = nSum + (Samples * vecCoeff2 * vecCoeff1);
    nDenom = nDenom + (vecCoeff2 * vecCoeff1); }}
    
    return nSum / nDenom;
}

float4 TexSharpenPass(float4 color)
{
    float3 calcSharpen = lumCoeff * A_SHARPEN_STRENGTH;

    float4 blurredColor = SampleBicubic(texcoord);
    float3 sharpenedColor = (color.rgb - blurredColor.rgb);

    float sharpenLuma = dot(sharpenedColor, calcSharpen);
    sharpenLuma = clamp(sharpenLuma, -B_SHARPEN_CLAMP, B_SHARPEN_CLAMP);

    color.rgb = color.rgb + sharpenLuma;
    color.a = AvgLuminance(color.rgb);

    return color;
}

/*------------------------------------------------------------------------------
                     [MAIN() & COMBINE PASS CODE SECTION]
------------------------------------------------------------------------------*/

void main()
{
    float4 color = Sample();

	color = LanczosScalerPass(color);
	//color = BicubicScalerPass(color);
	color = BilinearPass(color);
	color = TexSharpenPass(color);
	color = VibrancePass(color);

    SetOutput(color);
}
