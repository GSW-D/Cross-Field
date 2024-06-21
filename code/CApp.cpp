#include "CApp.h"

IMPLEMENT_APP(CApp);

bool CApp::OnInit()
{
    //ctor
    m_lpMainFrame = new CMainFrame(0L, _(_T("CrossField")));
    m_lpMainFrame->Show();
    return true;
}
