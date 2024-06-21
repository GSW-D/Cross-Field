#ifndef CGLCONTEXT_H
#define CGLCONTEXT_H

#include <wx/glcanvas.h>

class CGLContext : public wxGLContext
{
    public:
        CGLContext(wxGLCanvas *lpTarget);
    protected:
    private:
};

#endif // CGLCONTEXT_H
