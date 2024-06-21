#include "CGLContext.h"

CGLContext::CGLContext(wxGLCanvas *lpTarget):wxGLContext(lpTarget)
{
    //ctor
    SetCurrent(*lpTarget);
}
