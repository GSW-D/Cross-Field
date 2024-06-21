#include <wx/app.h>
#include "CMainFrame.h"
#ifndef CAPP_H
#define CAPP_H
class CApp : public wxApp
{
    public:
        virtual bool OnInit();
        CMainFrame *m_lpMainFrame;
    private:

};

#endif // CAPP_H
