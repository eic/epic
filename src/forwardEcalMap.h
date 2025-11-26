// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2025 Akio Ogawa
//==========================================================================
//  Implementation of forward calorimeter Mapping and numbering
//==========================================================================
//  Author: Akio Ogawa (BNL)
//==========================================================================

static const int mMaxNS       = 2;
static const int mMaxBlockId  = 1145;
static const int mMaxRowBlock = 39;
static const int mMaxColBlock = 19;

static const int mMaxTowerId  = 18320;
static const int mMaxRowTower = mMaxRowBlock * 4;
static const int mMaxColTower = mMaxColBlock * 4;

static const int mMaxFeebdId  = 594;
static const int mMaxRowFeebd = 39; //=mMaxRowBlock;
static const int mMaxColFeebd = 10; //=(mMaxColBlock+1)/2;

static const int mNColBlock[mMaxNS][mMaxRowBlock] = {
    {3,  6,  9,  11, 12, 13, 14, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 18, 18,
     18, 19, 19, 19, 18, 18, 18, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9,  6,  3},
    {3,  6,  9,  11, 12, 13, 14, 15, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 17, 17,
     17, 19, 19, 19, 18, 18, 18, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9,  6,  3}};

static const double mBackPlateZ = 362.0; //Global Z of backplate = start of Hcal
