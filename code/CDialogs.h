#include <wx/dialog.h>
#include <wx/button.h>
#include <wx/radiobut.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/sizer.h>
#include <wx/radiobox.h>
#include <wx/checkbox.h>
#include <sstream>
#include "CMeshInfoReserver.h"
//included by "CGLCanvas.h"
#ifndef CDIALOGS_H
#define CDIALOGS_H
namespace CDialogs
{
	class CTextureDialog : public wxDialog
	{
	public:

		CTextureDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("TextureOptions"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		bool m_bSelectedAreaTexture;
		bool m_bIncludingOtherVertices;
		bool m_bNormalizeTexture;
	protected:
		wxCheckBox* m_lpCheckBox1;
		wxCheckBox* m_lpCheckBox2;
		wxCheckBox* m_lpCheckBox3;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
		void OnChoosePart(wxCommandEvent &event);
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

    class CViewPortDialog : public wxDialog
    {
    public:
        typedef struct
        {
            double dblAxisX;
            double dblAxisY;
            double dblAxisZ;
            double dblRotateAngle;
            double dblDistance;
            double dblFarPlane;
            double dblNearPlane;
        } ViewPortData;

        CViewPortDialog( wxWindow* parent, ViewPortData *lpData, wxWindowID id = wxID_ANY, const wxString& title = wxT("ViewPortData"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( -1,-1 ), long style = wxDEFAULT_DIALOG_STYLE );
        ViewPortData *m_lpData;
	protected:
        wxStaticText* m_lpstaticTextRotationAxis;
        wxStaticText* m_lpStaticTextXLabel;
        wxTextCtrl* m_lpTextCtrlAxisX;
        wxStaticText* m_lpStaticTextYLabel;
        wxTextCtrl* m_lpTextCtrlAxisY;
        wxStaticText* m_lpStaticTextZLabel;
        wxTextCtrl* m_lpTextCtrlAxisZ;
        wxStaticText* m_lpStaticTextRotateAngle;
        wxTextCtrl* m_lpTextCtrlRotateAngle;
        wxStaticText* m_lpStaticTextDistance;
        wxTextCtrl* m_lpTextCtrlDistance;
        wxStaticText* m_StaticTextFarPlane;
        wxTextCtrl* m_lpTextCtrlFarPlane;
        wxStaticText* m_lpStaticTextNearPlane;
        wxTextCtrl* m_lpTextCtrlNearPlane;

        wxButton* m_lpButtonOK;
        wxButton* m_lpButton;

    private:
        void OnOK(wxCommandEvent &event);
        void OnCancel(wxCommandEvent &event);
        DECLARE_EVENT_TABLE()

    };

	class CAngleErrDialog : public wxDialog
	{
	public:

		CAngleErrDialog(wxWindow* lpParent, double *lpdblErr);
		~CAngleErrDialog();
		double m_lpdblErr[9];
		int m_iSelection;
	protected:
		wxRadioBox* m_lpRadioBox1;
		wxStaticText* m_staticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_staticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxStaticText* m_staticText3;
		wxTextCtrl* m_lpTextCtrl3;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		void OnRadioBox(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CLICDialog : public wxDialog
	{
	public:

		CLICDialog(wxWindow* parent);
		int iFieldChoice, iScaleOption;
		double dblTraceRadius, dblSingularityScale, dblMaxCrossGrad, dblMaxElectricField, dblMaxCurlField;
		bool bWithSingularities, bDiscrete;
	protected:
		wxRadioBox *m_lpRadioBox1, *m_lpRadioBox2;
		wxCheckBox *m_lpCheckBox1, *m_lpCheckBox2;
		wxStaticText *m_lpStaticText1, *m_lpStaticText2, *m_lpStaticText3;
		wxTextCtrl *m_lpTextCtrl1, *m_lpTextCtrl2, *m_lpTextCtrl3;
		wxButton *m_lpButtonOK, *m_lpButtonCancel;
	private:
		void OnRadioBox(wxCommandEvent &event);
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CFieldFileFormatChoiceDialog : public wxDialog
	{
	public:
		CFieldFileFormatChoiceDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = _("Choose Format"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		bool bWithIndex;
		bool bFrameField;
	protected:
		wxCheckBox* m_lpCheckBox1;
		wxCheckBox* m_lpCheckBox2;
		wxButton* m_lpOKButton;
		wxButton* m_lpCancelButton;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CChartInitDialog : public wxDialog
	{
	public:
		CChartInitDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Choose Initilize Mode"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
        int m_iChoice;
	protected:
		wxRadioBox* m_lpRadioBox;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
        void OnOK(wxCommandEvent &event);
        void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CEnergyDialog : public wxDialog
	{
	public:

		CEnergyDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		int m_iEnergyType;
	protected:
		wxRadioBox* m_lpRadioBox1;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CDirectionConstraintDialog : public wxDialog
	{
	public:

		CDirectionConstraintDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		double m_lpdblVector[3];
	protected:
		wxStaticText* m_lpStaticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_lpStaticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxStaticText* m_lpStaticText3;
		wxTextCtrl* m_lpTextCtrl3;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
        void OnOK(wxCommandEvent &event);
        void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	//class CEnergyDialog : public wxDialog
	//{
	//public:
	//	CEnergyDialog(wxWindow* parent);
	//	double dblRadius;
	//protected:
	//	wxStaticText* m_lpStaticText1;
	//	wxTextCtrl* m_lpTextCtrl1;
	//	wxStaticText* m_lpStaticText2;
	//	wxTextCtrl* m_lpTextCtrl2;
	//	wxStaticText* m_lpStaticText3;
	//	wxTextCtrl* m_lpTextCtrl3;
	//	wxButton* m_lpButtonRefresh;
	//	wxButton* m_lpButtonClose;
	//private:
	//	void OnOK(wxCommandEvent &event);
	//	void OnCancel(wxCommandEvent &event);
	//	DECLARE_EVENT_TABLE()
	//};
    class CIterDlg : public wxDialog
    {
    public:
        CIterDlg( wxWindow* parent, int nIterTimes = 0, wxWindowID id = wxID_ANY, const wxString& title = wxT("Iteration Options"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
        long m_nIterTimes;
	protected:
        wxStaticText* m_lpStaticTextIterTimesLabel;
        wxTextCtrl* m_lpTextCtrlIterTimesEdit;

        wxButton* m_lpButtonOK;
        wxButton* m_lpButtonCancel;
	private:
        void OnOK(wxCommandEvent &event);
        void OnCancel(wxCommandEvent &event);

        DECLARE_EVENT_TABLE()
    };

	class CFacetControlDialog : public wxDialog
	{

	public:
		long m_nNumFacets;
		bool m_bIsotropy;
		long m_nIterateTimes;
		CFacetControlDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Set Expected Facet Number"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
	protected:
		wxStaticText* m_lpStaticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_lpStaticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxCheckBox* m_lpCheckBox;
		wxButton* m_lpOKButton;
		wxButton* m_lpCancelButton;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()

	};

	class CAngleBiasDialog : public wxDialog
	{
	public:
		double m_dblAngleBias;
		CAngleBiasDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("InputAngleBias"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
	protected:
		wxStaticText* m_lpAngleBiasStaticText;
		wxTextCtrl* m_lpAngleBiasTextCtrl;
		wxButton* m_lpOKButton;
		wxButton* m_lpCancelButton;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CAdjustTypeDialog : public wxDialog
	{
	public:
		CAdjustTypeDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		int m_iLocalOptimizeType;
	protected:
		wxRadioBox* m_lpRadioBox;
		wxButton* m_lpOKButton;
		wxButton* m_lpCancelButton;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CGlobalAdjustDialog : public wxDialog
	{
	public:
		CGlobalAdjustDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("GlobalAdjustDialog"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		int m_iLocalOptimizeType;
		long m_nTestTimes;
		double m_dblChargeRadius;
	protected:
		wxRadioBox* m_lpRadioBox;
		wxStaticText* m_lpStaticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_lpStaticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxButton* m_lpOKButton;
		wxButton* m_lpCancelButton;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CEstimateDialog : public wxDialog
	{
	public:

		CEstimateDialog(wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxT("Estimate Parameters"), const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE);
		long m_nAspectRatioGroups, m_nAreaGroups;
		double m_dblAspectRatioUnit, m_dblLogAreaUnit;

	protected:
		wxStaticText* m_lpStaticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_lpStaticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxStaticText* m_lpStaticText3;
		wxTextCtrl* m_lpTextCtrl3;
		wxStaticText* m_lpStaticText4;
		wxTextCtrl* m_lpTextCtrl4;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;

	private:

		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};

	class CCountSingularityDialog : public wxDialog
	{
	public:
		CCountSingularityDialog(wxWindow* parent);
		double m_dblStartLevel, m_dblDecrement;
		long m_iTestTimes;
	protected:
		wxStaticText* m_lpStaticText1;
		wxTextCtrl* m_lpTextCtrl1;
		wxStaticText* m_lpStaticText2;
		wxTextCtrl* m_lpTextCtrl2;
		wxStaticText* m_lpStaticText3;
		wxTextCtrl* m_lpTextCtrl3;
		wxButton* m_lpButtonOK;
		wxButton* m_lpButtonCancel;
	private:
		void OnOK(wxCommandEvent &event);
		void OnCancel(wxCommandEvent &event);
		DECLARE_EVENT_TABLE()
	};
};



#endif // CDIALOGS_H
