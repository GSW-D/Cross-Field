#include "CQuadMeshEvaluator.h"
//included by "CMeshInfoReserver.h"


enum PROCESS_STATUS
{
	PS_EMPTY,
	PS_GEOMETRY_BLANK,
	PS_GEOMETRY_FILLED,
	PS_DIRECTION_SCATTER,
	PS_TRACK_PATH,
	PS_DIRECTION_DECOMPOSE,
	PS_DIRECTION_UNIFORM,
	PS_CREVASSE_GLOBAL,
	PS_PARAMETER_INT,
	PS_NEW_POINTS,
	PS_NEW_LINKAGE,
	PS_NEW_SCATTER,
	PS_NEW_UNIFORM,
	PS_NEW_COARSE,
	PS_NEW_OPTIMIZED
};

void m_fnStatusRetreat(PROCESS_STATUS psNewProcessStatus);

static PROCESS_STATUS& operator--(PROCESS_STATUS &psProcessStatus)
{
	int iValue = int(psProcessStatus);
	--iValue;
	psProcessStatus = PROCESS_STATUS(iValue);
	return psProcessStatus;
};
static PROCESS_STATUS& operator++(PROCESS_STATUS &psProcessStatus)
{
	int iValue = int(psProcessStatus);
	++iValue;
	psProcessStatus = PROCESS_STATUS(iValue);
	return psProcessStatus;
};

class CProcessStatusManager: public CMeshCenter
{
public:

	PROCESS_STATUS m_psProcessStatus;
	CProcessStatusManager();
	~CProcessStatusManager();
	void m_fnStatusRetreat(PROCESS_STATUS psNewProcessStatus);
private:
	void m_fnClearLinkage();
	void m_fnClearDensity();
	void m_fnScatterDirection();
	void m_fnClearCrevasse();
};