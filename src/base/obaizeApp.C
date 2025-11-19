#include "obaizeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
obaizeApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

obaizeApp::obaizeApp(const InputParameters & parameters) : MooseApp(parameters)
{
  obaizeApp::registerAll(_factory, _action_factory, _syntax);
}

obaizeApp::~obaizeApp() {}

void
obaizeApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<obaizeApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"obaizeApp"});
  Registry::registerActionsTo(af, {"obaizeApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
obaizeApp::registerApps()
{
  registerApp(obaizeApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
obaizeApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  obaizeApp::registerAll(f, af, s);
}
extern "C" void
obaizeApp__registerApps()
{
  obaizeApp::registerApps();
}
