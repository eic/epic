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

static const int mNTowerInBlock        = 4;     //a block contains 4x4 towers
static const double mBlockSize         = 10.0;  //size of block
static const double mBlockLength       = 17.0;  //size of block
static const double mSpaceBetweenBlock = 0.1;   //gap between blocks
static const double mBackPlateZ        = 362.0; //Global Z of backplate = start of Hcal
static const double mOffsetX[mMaxNS]   = {0.685 / 2.0,
                                          0.685 / 2.0}; //gap between north and south halves
static const double mOffsetY[mMaxNS]   = {
    0.0, 0.0}; //height offset between beamline and middle of detector
static const double mOffsetZ[mMaxNS] = {mBackPlateZ - mBlockLength, mBackPlateZ - mBlockLength};
//z position of front face of detector
static const double mOffsetXBeamPipe[mMaxNS] = {
    7.5, 20.05}; //3 rows at beamline height is shifted, also = insert width

