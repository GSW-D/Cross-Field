#include "CDialogs.h"
CDialogs::CTextureDialog::CTextureDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(4, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpCheckBox1 = new wxCheckBox(this, wxID_NOTOALL, wxT("Selected Area Texture"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpCheckBox1->SetValue(false);
	fgSizer1->Add(m_lpCheckBox1, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCheckBox2 = new wxCheckBox(this, wxID_ANY, wxT("Including Other Vertices"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpCheckBox2->Disable();
	fgSizer1->Add(m_lpCheckBox2, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCheckBox3 = new wxCheckBox(this, wxID_ANY, wxT("Normalize Texture"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer1->Add(m_lpCheckBox3, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CTextureDialog::OnChoosePart(wxCommandEvent &event)
{
	if (m_lpCheckBox1->GetValue())
	{
		m_lpCheckBox2->Enable();
	}
	else
	{
		m_lpCheckBox2->SetValue(true);
		m_lpCheckBox2->Disable();
	}
}

void CDialogs::CTextureDialog::OnOK(wxCommandEvent &event)
{
	m_bSelectedAreaTexture = m_lpCheckBox1->GetValue();
	m_bIncludingOtherVertices = m_lpCheckBox2->GetValue();
	m_bNormalizeTexture = m_lpCheckBox3->GetValue();
	EndModal(wxID_OK);
}

void CDialogs::CTextureDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CTextureDialog, wxDialog)
EVT_CHECKBOX(wxID_NOTOALL, CDialogs::CTextureDialog::OnChoosePart)
EVT_BUTTON(wxID_OK, CDialogs::CTextureDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CTextureDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CViewPortDialog::CViewPortDialog( wxWindow* parent, ViewPortData *lpData, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer( 7, 1, 0, 0 );

	m_lpstaticTextRotationAxis = new wxStaticText( this, wxID_ANY, wxT("RotationAxis"), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpstaticTextRotationAxis->Wrap( -1 );
	gSizer1->Add( m_lpstaticTextRotationAxis, 1, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL|wxALL, 5 );

	wxGridSizer* gSizer2;
	gSizer2 = new wxGridSizer( 1, 6, 0, 0 );

	m_lpStaticTextXLabel = new wxStaticText( this, wxID_ANY, wxT("x="), wxDefaultPosition, wxSize( -1,-1 ), 0 );
	m_lpStaticTextXLabel->Wrap( -1 );
	gSizer2->Add( m_lpStaticTextXLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxLEFT, 5 );

	m_lpTextCtrlAxisX = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( -1,-1 ), 0 );
	gSizer2->Add( m_lpTextCtrlAxisX, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpStaticTextYLabel = new wxStaticText( this, wxID_ANY, wxT("y="), wxDefaultPosition, wxSize( -1,-1 ), 0 );
	m_lpStaticTextYLabel->Wrap( -1 );
	gSizer2->Add( m_lpStaticTextYLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpTextCtrlAxisY = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer2->Add( m_lpTextCtrlAxisY, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpStaticTextZLabel = new wxStaticText( this, wxID_ANY, wxT("z="), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpStaticTextZLabel->Wrap( -1 );
	gSizer2->Add( m_lpStaticTextZLabel, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpTextCtrlAxisZ = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer2->Add( m_lpTextCtrlAxisZ, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	gSizer1->Add( gSizer2, 1, wxEXPAND|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxALL, 5 );

	wxGridSizer* gSizer3;
	gSizer3 = new wxGridSizer( 1, 2, 0, 0 );

	m_lpStaticTextRotateAngle = new wxStaticText( this, wxID_ANY, wxT("RotateAngle"), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpStaticTextRotateAngle->Wrap( -1 );
	gSizer3->Add( m_lpStaticTextRotateAngle, 0, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL|wxLEFT, 5 );

	m_lpTextCtrlRotateAngle = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer3->Add( m_lpTextCtrlRotateAngle, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	gSizer1->Add( gSizer3, 1, wxEXPAND|wxTOP|wxRIGHT|wxLEFT, 5 );

	wxGridSizer* gSizer5;
	gSizer5 = new wxGridSizer( 1, 2, 0, 0 );

	m_lpStaticTextDistance = new wxStaticText( this, wxID_ANY, wxT("Distance"), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpStaticTextDistance->Wrap( -1 );
	gSizer5->Add( m_lpStaticTextDistance, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxLEFT, 5 );

	m_lpTextCtrlDistance = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer5->Add( m_lpTextCtrlDistance, 0, wxALIGN_CENTER_VERTICAL|wxALIGN_CENTER_HORIZONTAL|wxRIGHT, 5 );

	gSizer1->Add( gSizer5, 1, wxEXPAND|wxRIGHT|wxLEFT, 5 );

	wxGridSizer* gSizer6;
	gSizer6 = new wxGridSizer( 1, 2, 0, 0 );

	m_StaticTextFarPlane = new wxStaticText( this, wxID_ANY, wxT("NearPlane"), wxDefaultPosition, wxDefaultSize, 0 );
	m_StaticTextFarPlane->Wrap( -1 );
	gSizer6->Add( m_StaticTextFarPlane, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxLEFT, 5 );

	m_lpTextCtrlFarPlane = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer6->Add( m_lpTextCtrlFarPlane, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	gSizer1->Add( gSizer6, 1, wxEXPAND|wxRIGHT|wxLEFT, 5 );

	wxGridSizer* gSizer7;
	gSizer7 = new wxGridSizer( 1, 2, 0, 0 );

	m_lpStaticTextNearPlane = new wxStaticText( this, wxID_ANY, wxT("FarPlane"), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpStaticTextNearPlane->Wrap( -1 );
	gSizer7->Add( m_lpStaticTextNearPlane, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxLEFT, 5 );

	m_lpTextCtrlNearPlane = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer7->Add( m_lpTextCtrlNearPlane, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxRIGHT, 5 );

	gSizer1->Add( gSizer7, 1, wxEXPAND|wxBOTTOM|wxRIGHT|wxLEFT, 5 );

	wxGridSizer* gSizer4;
	gSizer4 = new wxGridSizer( 1, 3, 0, 0 );


	gSizer4->Add( 0, 0, 1, wxEXPAND, 5 );

	m_lpButtonOK = new wxButton( this, wxID_OK, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer4->Add( m_lpButtonOK, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpButton = new wxButton( this, wxID_CANCEL, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer4->Add( m_lpButton, 0, wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL|wxBOTTOM|wxRIGHT|wxLEFT, 5 );

	gSizer1->Add( gSizer4, 1, wxEXPAND|wxALL, 5 );

	this->SetSizer( gSizer1 );
	this->Layout();
	gSizer1->Fit( this );

    m_lpData = lpData;
    if(m_lpData != NULL)
    {
        (*m_lpTextCtrlAxisX)<<(m_lpData->dblAxisX);
        (*m_lpTextCtrlAxisY)<<(m_lpData->dblAxisY);
        (*m_lpTextCtrlAxisZ)<<(m_lpData->dblAxisZ);
        (*m_lpTextCtrlRotateAngle)<<(m_lpData->dblRotateAngle);
        (*m_lpTextCtrlDistance)<<(m_lpData->dblDistance);
        (*m_lpTextCtrlFarPlane)<<(m_lpData->dblFarPlane);
        (*m_lpTextCtrlNearPlane)<<(m_lpData->dblNearPlane);
    }
}



void  CDialogs::CViewPortDialog::OnOK(wxCommandEvent &event)
{
    m_lpTextCtrlAxisX->GetLineText(0).ToDouble(&(m_lpData->dblAxisX));
    m_lpTextCtrlAxisY->GetLineText(0).ToDouble(&(m_lpData->dblAxisY));
    m_lpTextCtrlAxisZ->GetLineText(0).ToDouble(&(m_lpData->dblAxisZ));
    m_lpTextCtrlRotateAngle->GetLineText(0).ToDouble(&(m_lpData->dblRotateAngle));
    m_lpTextCtrlDistance->GetLineText(0).ToDouble(&(m_lpData->dblDistance));
    m_lpTextCtrlFarPlane->GetLineText(0).ToDouble(&(m_lpData->dblFarPlane));
    m_lpTextCtrlNearPlane->GetLineText(0).ToDouble(&(m_lpData->dblNearPlane));
    EndModal(wxID_OK);
}


void CDialogs::CViewPortDialog::OnCancel(wxCommandEvent &event)
{
    EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CViewPortDialog, wxDialog)
    EVT_BUTTON(wxID_OK, CDialogs::CViewPortDialog::OnOK)
    EVT_BUTTON(wxID_CANCEL, CDialogs::CViewPortDialog::OnCancel)
END_EVENT_TABLE()


CDialogs::CAngleErrDialog::CAngleErrDialog(wxWindow* lpParent, double *lpdblErr) : wxDialog(lpParent, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE)
{
	int i;
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(3, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxString m_radioBox1Choices[] = { wxT("U"), wxT("V"), wxT("Ave") };
	int m_radioBox1NChoices = sizeof(m_radioBox1Choices) / sizeof(wxString);
	m_lpRadioBox1 = new wxRadioBox(this, wxID_SELECTALL, wxT("wxRadioBox"), wxDefaultPosition, wxDefaultSize, m_radioBox1NChoices, m_radioBox1Choices, 3, wxRA_SPECIFY_COLS);
	m_lpRadioBox1->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox1, 0, wxALL, 5);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(3, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_staticText1 = new wxStaticText(this, wxID_ANY, wxT("min"), wxDefaultPosition, wxDefaultSize, 0);
	m_staticText1->Wrap(-1);
	fgSizer2->Add(m_staticText1, 0, wxALL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxString::FromDouble(lpdblErr[0]), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL, 5);

	m_staticText2 = new wxStaticText(this, wxID_ANY, wxT("ave"), wxDefaultPosition, wxDefaultSize, 0);
	m_staticText2->Wrap(-1);
	fgSizer2->Add(m_staticText2, 0, wxALL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxString::FromDouble(lpdblErr[1]), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL, 5);

	m_staticText3 = new wxStaticText(this, wxID_ANY, wxT("max"), wxDefaultPosition, wxDefaultSize, 0);
	m_staticText3->Wrap(-1);
	fgSizer2->Add(m_staticText3, 0, wxALL, 5);

	m_lpTextCtrl3 = new wxTextCtrl(this, wxID_ANY, wxString::FromDouble(lpdblErr[2]), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl3, 0, wxALL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
	for (i = 0; i < 9; ++i)
	{
		m_lpdblErr[i] = lpdblErr[i];
	}
}

CDialogs::CAngleErrDialog::~CAngleErrDialog()
{
}

void CDialogs::CAngleErrDialog::OnOK(wxCommandEvent &event)
{
	m_lpTextCtrl1->GetValue().ToDouble(&(m_lpdblErr[m_iSelection * 3]));
	m_lpTextCtrl3->GetValue().ToDouble(&(m_lpdblErr[m_iSelection * 3 + 2]));
	EndModal(wxID_OK);
}

void CDialogs::CAngleErrDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

void CDialogs::CAngleErrDialog::OnRadioBox(wxCommandEvent &event)
{
	m_iSelection = m_lpRadioBox1->GetSelection();
	m_lpTextCtrl1->SetLabelText(wxString::FromDouble(m_lpdblErr[m_iSelection * 3]));
	m_lpTextCtrl2->SetLabelText(wxString::FromDouble(m_lpdblErr[m_iSelection * 3 + 1]));
	m_lpTextCtrl3->SetLabelText(wxString::FromDouble(m_lpdblErr[m_iSelection * 3 + 2]));
	Refresh();
}

BEGIN_EVENT_TABLE(CDialogs::CAngleErrDialog, wxDialog)
	EVT_BUTTON(wxID_OK, CDialogs::CAngleErrDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CAngleErrDialog::OnCancel)
	EVT_RADIOBOX(wxID_SELECTALL, CDialogs::CAngleErrDialog::OnRadioBox)
END_EVENT_TABLE()

CDialogs::CLICDialog::CLICDialog(wxWindow* parent) : wxDialog(parent, wxID_ANY, wxT("LIC Options"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(7, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxString m_lpRadioBox1Choices[] = { wxT("Cross Field"), wxT("Cross Grad"), wxT("Electric Field"), wxT("Curl Field") };
	int m_lpRadioBox1NChoices = sizeof(m_lpRadioBox1Choices) / sizeof(wxString);
	m_lpRadioBox1 = new wxRadioBox(this, wxID_EDIT, wxT("Choose Field"), wxDefaultPosition, wxDefaultSize, m_lpRadioBox1NChoices, m_lpRadioBox1Choices, 3, wxRA_SPECIFY_COLS);
	m_lpRadioBox1->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(3, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("Trace Radius"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxT("10.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("Singularity Radius Scale"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxT("1.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText3 = new wxStaticText(this, wxID_ANY, wxT("Max Field Scale"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText3->Wrap(-1);
	fgSizer2->Add(m_lpStaticText3, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl3 = new wxTextCtrl(this, wxID_ANY, wxT("1.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl3, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxFlexGridSizer* fgSizer3;
	fgSizer3 = new wxFlexGridSizer(1, 2, 0, 0);
	fgSizer3->SetFlexibleDirection(wxBOTH);
	fgSizer3->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpCheckBox1 = new wxCheckBox(this, wxID_ANY, wxT("With Singularities"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpCheckBox1->SetValue(true);
	fgSizer3->Add(m_lpCheckBox1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL | wxEXPAND, 5);

	m_lpCheckBox2 = new wxCheckBox(this, wxID_ANY, wxT("Discrete"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer3->Add(m_lpCheckBox2, 0, wxALL, 5);

	fgSizer1->Add(fgSizer3, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxString m_lpRadioBox2Choices[] = { wxT("No Scale"), wxT("Scale Field"), wxT("Pure Scale") };
	int m_lpRadioBox2NChoices = sizeof(m_lpRadioBox2Choices) / sizeof(wxString);
	m_lpRadioBox2 = new wxRadioBox(this, wxID_ANY, wxT("Scale Options"), wxDefaultPosition, wxDefaultSize, m_lpRadioBox2NChoices, m_lpRadioBox2Choices, 3, wxRA_SPECIFY_COLS);
	m_lpRadioBox2->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox2, 0, wxALL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CLICDialog::OnRadioBox(wxCommandEvent &event)
{
	iFieldChoice = m_lpRadioBox1->GetSelection();
	if (iFieldChoice == 0)
	{
		m_lpTextCtrl3->SetLabelText(wxString::Format(_("%2.1lf"), 1.0));
	}
	if (iFieldChoice == 1)
	{
		m_lpTextCtrl3->SetLabelText(wxString::Format(_("%.3e"), dblMaxCrossGrad));
	}
	if (iFieldChoice == 2)
	{
		m_lpTextCtrl3->SetLabelText(wxString::Format(_("%.3e"), dblMaxElectricField));
	}
	if (iFieldChoice == 3)
	{
		m_lpTextCtrl3->SetLabelText(wxString::Format(_("%.3e"), dblMaxCurlField));
	}
}

void CDialogs::CLICDialog::OnOK(wxCommandEvent &event)
{
	iFieldChoice = m_lpRadioBox1->GetSelection();
	m_lpTextCtrl1->GetLineText(0).ToDouble(&dblTraceRadius);
	m_lpTextCtrl2->GetLineText(0).ToDouble(&dblSingularityScale);
	if (iFieldChoice == 1)
	{
		m_lpTextCtrl3->GetLineText(0).ToDouble(&dblMaxCrossGrad);
	}
	if (iFieldChoice == 2)
	{
		m_lpTextCtrl3->GetLineText(0).ToDouble(&dblMaxElectricField);
	}
	if (iFieldChoice == 3)
	{
		m_lpTextCtrl3->GetLineText(0).ToDouble(&dblMaxCurlField);
	}
	bWithSingularities = m_lpCheckBox1->GetValue();
	if (!bWithSingularities)
	{
		dblSingularityScale = 0.0;
	}
	bDiscrete = m_lpCheckBox2->GetValue();
	iScaleOption = m_lpRadioBox2->GetSelection();
	EndModal(wxID_OK);
}

void CDialogs::CLICDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CLICDialog, wxDialog)
	EVT_RADIOBOX(wxID_EDIT, CDialogs::CLICDialog::OnRadioBox)
	EVT_BUTTON(wxID_OK, CDialogs::CLICDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CLICDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CFieldFileFormatChoiceDialog::CFieldFileFormatChoiceDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(3, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpCheckBox1 = new wxCheckBox(this, wxID_ANY, wxT("With Index"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer1->Add(m_lpCheckBox1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCheckBox2 = new wxCheckBox(this, wxID_ANY, wxT("FrameField"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer1->Add(m_lpCheckBox2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpOKButton = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpOKButton, 0, wxALL, 5);

	m_lpCancelButton = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpCancelButton, 0, wxALL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CFieldFileFormatChoiceDialog::OnOK(wxCommandEvent &event)
{
	bWithIndex = m_lpCheckBox1->GetValue();
	bFrameField = m_lpCheckBox2->GetValue();
	EndModal(wxID_OK);
}

void CDialogs::CFieldFileFormatChoiceDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CFieldFileFormatChoiceDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CFieldFileFormatChoiceDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CFieldFileFormatChoiceDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CChartInitDialog::CChartInitDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxString m_lpRadioBoxChoices[] = { wxT("Extend"), wxT("Random") };
	int m_lpRadioBoxNChoices = sizeof(m_lpRadioBoxChoices) / sizeof(wxString);
	m_lpRadioBox = new wxRadioBox(this, wxID_ANY, wxT("InitializeMode"), wxDefaultPosition, wxDefaultSize, m_lpRadioBoxNChoices, m_lpRadioBoxChoices, 2, wxRA_SPECIFY_COLS);
	m_lpRadioBox->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}


void CDialogs::CChartInitDialog::OnOK(wxCommandEvent &event)
{
	m_iChoice = m_lpRadioBox->GetSelection();
	EndModal(wxID_OK);
}

void CDialogs::CChartInitDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}


BEGIN_EVENT_TABLE(CDialogs::CChartInitDialog, wxDialog)
	EVT_BUTTON(wxID_OK, CDialogs::CChartInitDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CChartInitDialog::OnCancel)
END_EVENT_TABLE()


CDialogs::CEnergyDialog::CEnergyDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);
	wxString m_lpRadioBox1Choices[] = { wxT("Smooth Energy"), wxT("Conformal Energy"), wxT("Curl Energy"), wxT("Discrete Energy"), wxT("MIQ Energy"), wxT("Align Energy") };
	int m_lpRadioBox1NChoices = sizeof(m_lpRadioBox1Choices) / sizeof(wxString);
	m_lpRadioBox1 = new wxRadioBox(this, wxID_ANY, wxT("Choose Energy Type"), wxDefaultPosition, wxDefaultSize, m_lpRadioBox1NChoices, m_lpRadioBox1Choices, 6, wxRA_SPECIFY_ROWS);
	m_lpRadioBox1->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox1, 0, wxALL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CEnergyDialog::OnOK(wxCommandEvent &event)
{
	m_iEnergyType = m_lpRadioBox1->GetSelection();
	EndModal(wxID_OK);
}

void CDialogs::CEnergyDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CEnergyDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CEnergyDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CEnergyDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CDirectionConstraintDialog::CDirectionConstraintDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(3, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("x"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxT("0.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("y"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxT("0.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL, 5);

	m_lpStaticText3 = new wxStaticText(this, wxID_ANY, wxT("z"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText3->Wrap(-1);
	fgSizer2->Add(m_lpStaticText3, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl3 = new wxTextCtrl(this, wxID_ANY, wxT("1.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl3, 0, wxALL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CDirectionConstraintDialog::OnOK(wxCommandEvent &event)
{
	m_lpTextCtrl1->GetLineText(0).ToDouble(m_lpdblVector);
	m_lpTextCtrl2->GetLineText(0).ToDouble(m_lpdblVector + 1);
	m_lpTextCtrl3->GetLineText(0).ToDouble(m_lpdblVector + 2);
	EndModal(wxID_OK);
}

void CDialogs::CDirectionConstraintDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CDirectionConstraintDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CDirectionConstraintDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CDirectionConstraintDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CIterDlg::CIterDlg( wxWindow* parent, int nIterTimes, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxDialog( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer( 2, 1, 0, 0 );

	wxGridSizer* gSizer2;
	gSizer2 = new wxGridSizer( 1, 2, 0, 0 );

	m_lpStaticTextIterTimesLabel = new wxStaticText( this, wxID_ANY, wxT("Iteration Times"), wxDefaultPosition, wxDefaultSize, 0 );
	m_lpStaticTextIterTimesLabel->Wrap( -1 );
	gSizer2->Add( m_lpStaticTextIterTimesLabel, 0, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpTextCtrlIterTimesEdit = new wxTextCtrl( this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0 );
	gSizer2->Add( m_lpTextCtrlIterTimesEdit, 0, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	gSizer1->Add( gSizer2, 1, wxEXPAND, 5 );

	wxGridSizer* gSizer3;
	gSizer3 = new wxGridSizer( 1, 3, 0, 0 );


	gSizer3->Add( 0, 0, 1, wxEXPAND, 5 );

	m_lpButtonOK = new wxButton( this, wxID_OK, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer3->Add( m_lpButtonOK, 0, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	m_lpButtonCancel = new wxButton( this, wxID_CANCEL, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0 );
	gSizer3->Add( m_lpButtonCancel, 0, wxALL|wxALIGN_CENTER_HORIZONTAL|wxALIGN_CENTER_VERTICAL, 5 );

	gSizer1->Add( gSizer3, 1, wxEXPAND, 5 );

	this->SetSizer( gSizer1 );
	this->Layout();
	gSizer1->Fit( this );
    (*m_lpTextCtrlIterTimesEdit) << nIterTimes;
}



void CDialogs::CIterDlg::OnOK(wxCommandEvent &event)
{
    m_lpTextCtrlIterTimesEdit->GetLineText(0).ToLong(&m_nIterTimes);
    EndModal(wxID_OK);
}

void CDialogs::CIterDlg::OnCancel(wxCommandEvent &event)
{
    EndModal(wxID_CANCEL);
}


BEGIN_EVENT_TABLE(CDialogs::CIterDlg, wxDialog)
    EVT_BUTTON(wxID_OK, CDialogs::CIterDlg::OnOK)
    EVT_BUTTON(wxID_CANCEL, CDialogs::CIterDlg::OnCancel)
END_EVENT_TABLE()

CDialogs::CFacetControlDialog::CFacetControlDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(3, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(2, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("Expected Facet Number"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxT("5000"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("Iterate Times"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxT("1"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	fgSizer1->Add(fgSizer2, 1, wxEXPAND, 5);

	m_lpCheckBox = new wxCheckBox(this, wxID_ANY, wxT("Isotropy"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpCheckBox->SetValue(true);
	fgSizer1->Add(m_lpCheckBox, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpOKButton = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpOKButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCancelButton = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpCancelButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);

	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CFacetControlDialog::OnOK(wxCommandEvent &event)
{
	m_lpTextCtrl1->GetLineText(0).ToLong(&m_nNumFacets);
	m_lpTextCtrl2->GetLineText(0).ToLong(&m_nIterateTimes);
	m_bIsotropy = m_lpCheckBox->GetValue();
	EndModal(wxID_OK);
}

void CDialogs::CFacetControlDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CFacetControlDialog, wxDialog)
	EVT_BUTTON(wxID_OK, CDialogs::CFacetControlDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CFacetControlDialog::OnCancel)
END_EVENT_TABLE()


CDialogs::CAngleBiasDialog::CAngleBiasDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(2, 2, 0, 0);

	m_lpAngleBiasStaticText = new wxStaticText(this, wxID_ANY, wxT("AngleBias"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpAngleBiasStaticText->Wrap(-1);
	gSizer1->Add(m_lpAngleBiasStaticText, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpAngleBiasTextCtrl = new wxTextCtrl(this, wxID_ANY, wxEmptyString, wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpAngleBiasTextCtrl, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpOKButton = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpOKButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCancelButton = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpCancelButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(gSizer1);
	this->Layout();
	gSizer1->Fit(this);

	this->Centre(wxBOTH);
}




void CDialogs::CAngleBiasDialog::OnOK(wxCommandEvent &event)
{
	m_lpAngleBiasTextCtrl->GetLineText(0).ToDouble(&m_dblAngleBias);
	EndModal(wxID_OK);
}

void CDialogs::CAngleBiasDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CAngleBiasDialog, wxDialog)
	EVT_BUTTON(wxID_OK, CDialogs::CAngleBiasDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CAngleBiasDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CAdjustTypeDialog::CAdjustTypeDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxString m_lpRadioBoxChoices[] = { wxT("Chart Direction"), wxT("Conformal Coefficient") };
	int m_lpRadioBoxNChoices = sizeof(m_lpRadioBoxChoices) / sizeof(wxString);
	m_lpRadioBox = new wxRadioBox(this, wxID_ANY, wxT("Field Selection"), wxDefaultPosition, wxDefaultSize, m_lpRadioBoxNChoices, m_lpRadioBoxChoices, 2, wxRA_SPECIFY_COLS);
	m_lpRadioBox->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox, 0, wxALL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpOKButton = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpOKButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpCancelButton = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpCancelButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CAdjustTypeDialog::OnOK(wxCommandEvent &event)
{
	m_iLocalOptimizeType = m_lpRadioBox->GetSelection();
	EndModal(wxID_OK);
}

void CDialogs::CAdjustTypeDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CAdjustTypeDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CAdjustTypeDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CAdjustTypeDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CGlobalAdjustDialog::CGlobalAdjustDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(3, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxString m_lpRadioBoxChoices[] = { wxT("Chart Direction"), wxT("Conformal Coefficient") };
	int m_lpRadioBoxNChoices = sizeof(m_lpRadioBoxChoices) / sizeof(wxString);
	m_lpRadioBox = new wxRadioBox(this, wxID_ANY, wxT("Field Selection"), wxDefaultPosition, wxDefaultSize, m_lpRadioBoxNChoices, m_lpRadioBoxChoices, 2, wxRA_SPECIFY_COLS);
	m_lpRadioBox->SetSelection(0);
	fgSizer1->Add(m_lpRadioBox, 0, wxALL, 5);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(2, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("Test Times"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxT("3"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("Charge Radius"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxT("-2"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer3;
	gSizer3 = new wxGridSizer(1, 2, 0, 0);

	m_lpOKButton = new wxButton(this, wxID_OK, wxT("OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer3->Add(m_lpOKButton, 0, wxALL, 5);

	m_lpCancelButton = new wxButton(this, wxID_CANCEL, wxT("Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer3->Add(m_lpCancelButton, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer3, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}



void CDialogs::CGlobalAdjustDialog::OnOK(wxCommandEvent &event)
{
	m_iLocalOptimizeType = m_lpRadioBox->GetSelection();
	m_lpTextCtrl1->GetLineText(0).ToLong(&m_nTestTimes);
	m_lpTextCtrl2->GetLineText(0).ToDouble(&m_dblChargeRadius);
	EndModal(wxID_OK);
}

void CDialogs::CGlobalAdjustDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CGlobalAdjustDialog, wxDialog)
	EVT_BUTTON(wxID_OK, CDialogs::CGlobalAdjustDialog::OnOK)
	EVT_BUTTON(wxID_CANCEL, CDialogs::CGlobalAdjustDialog::OnCancel)
END_EVENT_TABLE()


CDialogs::CEstimateDialog::CEstimateDialog(wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style) : wxDialog(parent, id, title, pos, size, style)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(4, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("Aspect Ratio Unit"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, _("0.125"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("Aspect Ratio Groups"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, _("7"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText3 = new wxStaticText(this, wxID_ANY, wxT("Area Unit"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText3->Wrap(-1);
	fgSizer2->Add(m_lpStaticText3, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl3 = new wxTextCtrl(this, wxID_ANY, _("0.25"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl3, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpStaticText4 = new wxStaticText(this, wxID_ANY, wxT("Area Groups"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText4->Wrap(-1);
	fgSizer2->Add(m_lpStaticText4, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpTextCtrl4 = new wxTextCtrl(this, wxID_ANY, _("7"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl4, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND | wxALIGN_CENTER_HORIZONTAL | wxALIGN_CENTER_VERTICAL, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CEstimateDialog::OnOK(wxCommandEvent &event)
{
	m_lpTextCtrl1->GetLineText(0).ToDouble(&m_dblAspectRatioUnit);
	m_lpTextCtrl2->GetLineText(0).ToLong(&m_nAspectRatioGroups);
	m_lpTextCtrl3->GetLineText(0).ToDouble(&m_dblLogAreaUnit);
	m_lpTextCtrl4->GetLineText(0).ToLong(&m_nAreaGroups);
	EndModal(wxID_OK);
}

void CDialogs::CEstimateDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CEstimateDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CEstimateDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CEstimateDialog::OnCancel)
END_EVENT_TABLE()

CDialogs::CCountSingularityDialog::CCountSingularityDialog(wxWindow* parent) : wxDialog(parent, wxID_ANY, wxT("Count Singularities"), wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE)
{
	this->SetSizeHints(wxDefaultSize, wxDefaultSize);

	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer(2, 1, 0, 0);
	fgSizer1->SetFlexibleDirection(wxBOTH);
	fgSizer1->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	wxFlexGridSizer* fgSizer2;
	fgSizer2 = new wxFlexGridSizer(3, 2, 0, 0);
	fgSizer2->SetFlexibleDirection(wxBOTH);
	fgSizer2->SetNonFlexibleGrowMode(wxFLEX_GROWMODE_SPECIFIED);

	m_lpStaticText1 = new wxStaticText(this, wxID_ANY, wxT("Start Level"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText1->Wrap(-1);
	fgSizer2->Add(m_lpStaticText1, 0, wxALL, 5);

	m_lpTextCtrl1 = new wxTextCtrl(this, wxID_ANY, wxT("-1.0"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl1, 0, wxALL, 5);

	m_lpStaticText2 = new wxStaticText(this, wxID_ANY, wxT("Decrement"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText2->Wrap(-1);
	fgSizer2->Add(m_lpStaticText2, 0, wxALL, 5);

	m_lpTextCtrl2 = new wxTextCtrl(this, wxID_ANY, wxT("0.05"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl2, 0, wxALL, 5);

	m_lpStaticText3 = new wxStaticText(this, wxID_ANY, wxT("Test Times"), wxDefaultPosition, wxDefaultSize, 0);
	m_lpStaticText3->Wrap(-1);
	fgSizer2->Add(m_lpStaticText3, 0, wxALL, 5);

	m_lpTextCtrl3 = new wxTextCtrl(this, wxID_ANY, wxT("60"), wxDefaultPosition, wxDefaultSize, 0);
	fgSizer2->Add(m_lpTextCtrl3, 0, wxALL, 5);


	fgSizer1->Add(fgSizer2, 1, wxEXPAND, 5);

	wxGridSizer* gSizer1;
	gSizer1 = new wxGridSizer(1, 2, 0, 0);

	m_lpButtonOK = new wxButton(this, wxID_OK, wxT("&OK"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonOK, 0, wxALL, 5);

	m_lpButtonCancel = new wxButton(this, wxID_CANCEL, wxT("&Cancel"), wxDefaultPosition, wxDefaultSize, 0);
	gSizer1->Add(m_lpButtonCancel, 0, wxALL, 5);


	fgSizer1->Add(gSizer1, 1, wxEXPAND, 5);


	this->SetSizer(fgSizer1);
	this->Layout();
	fgSizer1->Fit(this);

	this->Centre(wxBOTH);
}

void CDialogs::CCountSingularityDialog::OnOK(wxCommandEvent &event)
{
	m_lpTextCtrl1->GetLineText(0).ToDouble(&m_dblStartLevel);
	m_lpTextCtrl2->GetLineText(0).ToDouble(&m_dblDecrement);
	m_lpTextCtrl3->GetLineText(0).ToLong(&m_iTestTimes);
	EndModal(wxID_OK);
}

void CDialogs::CCountSingularityDialog::OnCancel(wxCommandEvent &event)
{
	EndModal(wxID_CANCEL);
}

BEGIN_EVENT_TABLE(CDialogs::CCountSingularityDialog, wxDialog)
EVT_BUTTON(wxID_OK, CDialogs::CCountSingularityDialog::OnOK)
EVT_BUTTON(wxID_CANCEL, CDialogs::CCountSingularityDialog::OnCancel)
END_EVENT_TABLE()
