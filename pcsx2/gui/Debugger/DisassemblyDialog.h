/*  PCSX2 - PS2 Emulator for PCs
 *  Copyright (C) 2002-2014  PCSX2 Dev Team
 *
 *  PCSX2 is free software: you can redistribute it and/or modify it under the terms
 *  of the GNU Lesser General Public License as published by the Free Software Found-
 *  ation, either version 3 of the License, or (at your option) any later version.
 *
 *  PCSX2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along with PCSX2.
 *  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#include <wx/wx.h>
#include <wx/notebook.h>
#include "App.h"

#include "CtrlDisassemblyView.h"
#include "CtrlRegisterList.h"
#include "CtrlMemView.h"
#include "DebugEvents.h"
#include "DebuggerLists.h"

class DebuggerHelpDialog: public wxDialog
{
public:
	DebuggerHelpDialog(wxWindow* parent);
};

class CpuTabPage: public wxPanel
{
public:
	CpuTabPage(wxWindow* parent, DebugInterface* _cpu);
	DebugInterface* getCpu() { return cpu; };
	CtrlDisassemblyView* getDisassembly() { return disassembly; };
	CtrlRegisterList* getRegisterList() { return registerList; };
	CtrlMemView* getMemoryView() { return memory; };
	wxNotebook* getBottomTabs() { return bottomTabs; };
	void update();
	void showMemoryView() { setBottomTabPage(memory); };
	void loadCycles();
	void reloadSymbolMap();
	u32 getStepOutAddress();

	void listBoxHandler(wxCommandEvent& event);
	DECLARE_EVENT_TABLE()
private:
	void setBottomTabPage(wxWindow* win);
	void postEvent(wxEventType type, int value);

	DebugInterface* cpu;
	CtrlDisassemblyView* disassembly;
	CtrlRegisterList* registerList;
	wxListBox* functionList;
	CtrlMemView* memory;
	wxNotebook* bottomTabs;
	wxNotebook* leftTabs;
	BreakpointList* breakpointList;
	wxStaticText* cyclesText;
	ThreadList* threadList;
	StackFramesList* stackFrames;
	u32 lastCycles;
};

class DisassemblyDialog : public wxFrame
{
public:
	DisassemblyDialog( wxWindow* parent=NULL );
	virtual ~DisassemblyDialog() throw() {}

	static wxString GetNameStatic() { return L"DisassemblyDialog"; }
	wxString GetDialogName() const { return GetNameStatic(); }
	
	void update();
	void reset();
	void setDebugMode(bool debugMode, bool switchPC);
	
#ifdef WIN32
	WXLRESULT MSWWindowProc(WXUINT nMsg, WXWPARAM wParam, WXLPARAM lParam);
#endif

	DECLARE_EVENT_TABLE()
protected:
	void onBreakRunClicked(wxCommandEvent& evt);
	void onStepOverClicked(wxCommandEvent& evt);
	void onStepIntoClicked(wxCommandEvent& evt);
	void onStepOutClicked(wxCommandEvent& evt);
	void onDebuggerEvent(wxCommandEvent& evt);
	void onPageChanging(wxCommandEvent& evt);
	void onBreakpointClick(wxCommandEvent& evt);
	void onClose(wxCloseEvent& evt);
	void stepOver();
	void stepInto();
	void stepOut();
	void gotoPc();
private:
	CpuTabPage* eeTab;
	CpuTabPage* iopTab;
	CpuTabPage* currentCpu;

	wxBoxSizer* topSizer;
	wxStatusBar* statusBar;
	wxButton* breakRunButton;
	wxButton* stepIntoButton;
	wxButton* stepOverButton;
	wxButton* stepOutButton;
	wxButton* breakpointButton;
};