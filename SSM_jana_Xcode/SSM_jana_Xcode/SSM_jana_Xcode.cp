/*
 *  SSM_jana_Xcode.cp
 *  SSM_jana_Xcode
 *
 *  Created by Lipkova on 03/04/17.
 *  Copyright (c) 2017 Lipkova. All rights reserved.
 *
 */

#include "SSM_jana_Xcode.h"
#include "SSM_jana_XcodePriv.h"

CFStringRef SSM_jana_XcodeUUID(void)
{
	CSSM_jana_Xcode* theObj = new CSSM_jana_Xcode;
	return theObj->UUID();
}

CFStringRef CSSM_jana_Xcode::UUID()
{
	return CFSTR("0001020304050607");
}
