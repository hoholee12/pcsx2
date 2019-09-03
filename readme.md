PCSX2 1.3.1 BGI(Better GameIndex) Version
fork of https://github.com/PCSX2/pcsx2


things added to the gameindex.dbf
=================================
FramerateNTSC, FrameratePAL = FixedInt100

VUCycleSteal = u8

EECycleRate = u8

vuThread = bool

VsyncEnable = bool

fastCDVD = bool

also in the source
==================
pcsx2iidxtuner https://github.com/hoholee12/pcsx2iidxtuner read howto.txt

base: https://github.com/PCSX2/pcsx2/tree/7b1214849acba37fd66693359d1e3fced001ed78
i chose this build because its the most stable and fastest version ive ever used.
cdvd portion is updated to 1.4.0 (cso, gzip support)
