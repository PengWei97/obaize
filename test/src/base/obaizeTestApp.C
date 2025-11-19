//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "obaizeTestApp.h"
#include "obaizeApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
obaizeTestApp::validParams()
{
  InputParameters params = obaizeApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

obaizeTestApp::obaizeTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  obaizeTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

obaizeTestApp::~obaizeTestApp() {}

void
obaizeTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  obaizeApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"obaizeTestApp"});
    Registry::registerActionsTo(af, {"obaizeTestApp"});
  }
}

void
obaizeTestApp::registerApps()
{
  registerApp(obaizeApp);
  registerApp(obaizeTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
obaizeTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  obaizeTestApp::registerAll(f, af, s);
}
extern "C" void
obaizeTestApp__registerApps()
{
  obaizeTestApp::registerApps();
}
