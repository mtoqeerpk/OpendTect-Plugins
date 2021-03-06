#include "uigeopackageexportmod.h"

#include "uimenu.h"
#include "uiodmenumgr.h"
#include "uitoolbar.h"
#include "uistring.h"
#include "odplugin.h"
#include "survinfo.h"
#include "uimsg.h"
#include "coordsystem.h"

#include "uigeopackageexportmainwin.h"


mDefODPluginInfo(uiGeopackageExport)
{
    mDefineStaticLocalObject( PluginInfo, retpi,(
    "Geopackage Export plugin",
    "",
    "Wayne Mogg",
    "1.0",
    "Export OpendTect survey information to a Geopackage database."));
    return &retpi;
}

class uiGeopackageExportMgr :  public CallBacker
{ mODTextTranslationClass(uiGeopackageExportMgr)
public:
			uiGeopackageExportMgr(uiODMain*);
			~uiGeopackageExportMgr();

    uiODMain*		appl_;
    uiGeopackageExportMainWin*	dlg_;

    void		updateToolBar(CallBacker*);
    void		updateMenu(CallBacker*);
    void		exportDialog(CallBacker*);
};


uiGeopackageExportMgr::uiGeopackageExportMgr( uiODMain* a )
	: appl_(a)
	, dlg_(0)
{
    mAttachCB( appl_->menuMgr().dTectTBChanged, uiGeopackageExportMgr::updateToolBar );
    mAttachCB( appl_->menuMgr().dTectMnuChanged, uiGeopackageExportMgr::updateMenu );
    updateToolBar(0);
    updateMenu(0);
}


uiGeopackageExportMgr::~uiGeopackageExportMgr()
{
    detachAllNotifiers();
}


void uiGeopackageExportMgr::updateToolBar( CallBacker* )
{
}


void uiGeopackageExportMgr::updateMenu( CallBacker* )
{
    delete dlg_;
    dlg_ = 0;
    uiAction* newitem = new uiAction( tr("Geopackage Export"), mCB(this,uiGeopackageExportMgr,exportDialog));
    appl_->menuMgr().getBaseMnu( uiODApplMgr::Exp )->insertItem( newitem );
}


void uiGeopackageExportMgr::exportDialog( CallBacker* )
{

    if ( !dlg_ ) {
        dlg_ = new uiGeopackageExportMainWin( appl_ );
        if (!dlg_->checkCRSdefined()) {
            SurveyInfo* si = const_cast<SurveyInfo*>( &SI() );
            uiString msg( tr("Geopackage Export requires a projected Coordinate Reference System (CRS)."
            " The survey '%1' currently has a CRS set to: '%2' which is not a projected CRS."
            " You can set the survey CRS using the Survey-Select/Setup menu item.").arg(si->name()).arg(si->getCoordSystem()->factoryKeyword()) );
            uiMSG().message(msg);
            dlg_->close();
            return;
        }
    }

    dlg_->show();
    dlg_->raise();
}


mDefODInitPlugin(uiGeopackageExport)
{
    mDefineStaticLocalObject( PtrMan<uiGeopackageExportMgr>, theinst_, = 0 );
    if ( theinst_ ) return 0;

    theinst_ = new uiGeopackageExportMgr( ODMainWin() );
    if ( !theinst_ )
        return "Cannot instantiate Geopackage Export plugin";
    
     return 0;
}
