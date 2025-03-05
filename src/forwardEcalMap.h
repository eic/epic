static const int mMaxNS=2;

static const int mMaxBlockId=1145;
static const int mMaxRowBlock=39;
static const int mMaxColBlock=19;

static const int mMaxTowerId=18320;
static const int mMaxRowTower=mMaxRowBlock*4;
static const int mMaxColTower=mMaxColBlock*4;

static const int mMaxFeebdId=594;
static const int mMaxRowFeebd=39; //=mMaxRowBlock;  
static const int mMaxColFeebd=10; //=(mMaxColBlock+1)/2; 
  
static const int mNColBlock[mMaxNS][mMaxRowBlock]={
  { 3, 6, 9,11,12,13,14,15,16,16,17,17,18,18,18,19,19,19,18,18,18,19,19,19,18,18,18,17,17,16,16,15,14,13,12,11, 9, 6, 3},
  { 3, 6, 9,11,12,13,14,15,16,16,17,17,18,18,18,19,19,19,17,17,17,19,19,19,18,18,18,17,17,16,16,15,14,13,12,11, 9, 6, 3}
};

static const int    mNTowerInBlock=4;                   //a block contains 4x4 towers
static const double mBlockSize=10.0;                    //size of block 
static const double mSpaceBetweenBlock=0.0254;          //gap between blocks
static const double mOffsetX[mMaxNS]={0.1,0.1};         //gap between north and south halves
static const double mOffsetY[mMaxNS]={0.0,0.0};         //height offset between beamline and middle of detector  
static const double mOffsetZ[mMaxNS]={300.0,300.0};     //z position of front face of detector
static const double mOffsetXBeamPipe[mMaxNS]={7.5,20.0};//3 rows at beamline height is shifted

class forwardEcalMap{
private:
  int mTowerNS [mMaxTowerId];
  int mTowerRow[mMaxTowerId];
  int mTowerCol[mMaxTowerId];
  int mTowerId [mMaxNS][mMaxRowTower][mMaxColTower];
  int mBlockNS [mMaxBlockId];
  int mBlockRow[mMaxBlockId];
  int mBlockCol[mMaxBlockId];
  int mBlockId [mMaxNS][mMaxRowBlock][mMaxColBlock];
  int mFeebdNS [mMaxFeebdId];
  int mFeebdRow[mMaxFeebdId];
  int mFeebdCol[mMaxFeebdId];
  int mFeebdId [mMaxNS][mMaxRowFeebd][mMaxColFeebd];

public:
  forwardEcalMap() {createMap();}
  ~forwardEcalMap() {};
  
  int maxNS(){return mMaxNS;}
  int maxRowBlock(){return mMaxRowBlock;}
  int maxColBlock(){return mMaxColBlock;}
  int nColBlock(int ns, int row){return mNColBlock[ns][row];}
  double blockSize(){return mBlockSize;}
  double spaceBetweenBlock(){return mSpaceBetweenBlock;}
  double offsetX(int ns){return mOffsetX[ns];}
  double offsetY(int ns){return mOffsetY[ns];}
  double offsetZ(int ns){return mOffsetZ[ns];}
  int blockNS (int bid) {return mBlockNS [bid];}
  int blockRow(int bid) {return mBlockRow[bid];}
  int blockCol(int bid) {return mBlockCol[bid];}
  int blockId(int ns, int rowB, int colB) {return mBlockId[ns][rowB][colB];}

  // Block of 4x4 towers is the basic unit
  double xBlock(int ns, int row, int col){
    if(row>=18 && row<=20){
      return (1-2*ns)*(mOffsetX[ns]+mOffsetXBeamPipe[ns]+(mBlockSize+mSpaceBetweenBlock)*(col+0.5));
    }else{
      return (1-2*ns)*(mOffsetX[ns]+(mBlockSize+mSpaceBetweenBlock)*(col+0.5));
    }
  }
  double yBlock(int ns, int row){
    return mOffsetY[ns]+(mBlockSize+mSpaceBetweenBlock)*(mMaxRowBlock/2.0-row-0.5);
  }
  double doesBlockExist(int ns, int row, int col){
    if(0<=row && row<mMaxRowBlock && 0<=col && col<mNColBlock[ns][row]) return 1;
    return 0;
  }
  
  // Each block contains 4x4 tower  
  int maxRowTower(){return mMaxRowBlock*mNTowerInBlock;}
  int maxColTower(){return mMaxColBlock*mNTowerInBlock;}
  int nColTower(int ns, int row){return mNColBlock[ns][rowBlock(row)]*mNTowerInBlock;}
  double towerSize(){return mBlockSize/double(mNTowerInBlock);}
  int colBlock(int colTower){return colTower/mNTowerInBlock;}
  int rowBlock(int rowTower){return rowTower/mNTowerInBlock;}
  double xTower(int ns, int row, int col){return  xBlock(ns, rowBlock(row), colBlock(col)) + (1-2*ns)*(col%mNTowerInBlock - 1.5)*towerSize();}
  double yTower(int ns, int row){return  yBlock(ns, rowBlock(row)) + (1-2*ns)*(row%mNTowerInBlock - 1.5)*towerSize();}
  double doesTowerExist(int ns, int row, int col){return doesBlockExist(ns,rowBlock(row),colBlock(col));}
  int towerNS (int tid) {return mTowerNS [tid];}
  int towerRow(int tid) {return mTowerRow[tid];}
  int towerCol(int tid) {return mTowerCol[tid];}
  int towerId(int ns, int rowT, int colT) {return mTowerId[ns][rowT][colT];}

  // Each FEB contains 1x2 blocks
  int maxRowFeebd(){return mMaxRowBlock;}
  int maxColFeebd(){return (mMaxColBlock+1)/2;}
  int nColFeebd(int ns, int row){return (mNColBlock[ns][row]+1)/2;}
  double feebdSizeX(){return mBlockSize*2;}
  double feebdSizeY(){return mBlockSize;}
  int colFeebd(int colBlock){return colBlock/2;}
  int rowFeebd(int rowBlock){return rowBlock;}
  double xFeebd(int ns, int row, int col){return  xBlock(ns, row, col*2) + (0.5-ns)*blockSize();}
  double yFeebd(int ns, int row){return  yBlock(ns, row);}
  double doesFeebdExist(int ns, int row, int col){return doesBlockExist(ns,row,col*2);}
  int feebdNS (int fid) {return mFeebdNS [fid];}
  int feebdRow(int fid) {return mFeebdRow[fid];}
  int feebdCol(int fid) {return mFeebdCol[fid];}
  int feebdId(int ns, int rowF, int colF) {return mFeebdId[ns][rowF][colF];}
  
  // Power Group
  static const int mNPowerGroup=34;
  int powerGroup(int rowB, int colB){
    if     (rowB== 0) {return 1;}
    else if(rowB== 1) {return 1;}
    else if(rowB== 2) {if(colB<2) {return 2;} else {return 1;} }
    else if(rowB== 3) {return 2;}
    else if(rowB== 4) {if(colB<2) {return 4;} else if(colB<8){return 3;} else {return 2;}}
    else if(rowB== 5) {if(colB<2) {return 4;} else {return 3;}}
    else if(rowB>= 6  && rowB<=14) {return rowB-2;}
    else if(rowB>=15 && rowB<=17) {if(colB<18) {return rowB-2;} else {return rowB-10;} }
    else if(rowB>=18 && rowB<=20) {return rowB-2;}
    else if(rowB>=21 && rowB<=23) {if(colB<18) {return rowB-2;} else {return rowB+6;} }
    else if(rowB>=24 && rowB<=32) {return rowB-2;}
    else if(rowB==33) {if(colB<2) {return 30;} else {return 31;}}
    else if(rowB==34) {if(colB<2) {return 30;} else if(colB<8){return 31;} else {return 32;}}
    else if(rowB==35) {return 32;}
    else if(rowB==36) {if(colB<2) {return 32;} else {return 33;} }
    else if(rowB==37) {return 33;}
    else if(rowB==38) {return 33;}
    return -1;
  }  
  int powerGroupColor(int pg){
    int pgc[mNPowerGroup]={
      kOrange,kGreen-9,kViolet-9,kCyan-9,
      kGreen+2,kAzure+7,kViolet,kMagenta-9,kGreen-9,
      kViolet-9,kCyan-9,kGray,
      kOrange,kMagenta-9,kGreen-9,
      kViolet-9,kCyan-9,kViolet-9,
      kGreen-9,kMagenta-9,kOrange,
      kGray,kCyan-9,kViolet-9,
      kGreen-9,kMagenta-9,kViolet,kAzure+7,kGreen+2,
      kCyan-9,kViolet-9,kGreen-9,kOrange};
    return pgc[pg-1];
  }
  
  void createMap(){
    int tid=0, bid=0, fid=0;
    for(int ns=0; ns<maxNS(); ns++){
      for(int row=0; row<maxRowTower(); row++){
	for(int col=0; col<maxColTower(); col++){
	  if(doesTowerExist(ns,row,col)){
	    mTowerNS [tid]=ns;
	    mTowerRow[tid]=row;
	    mTowerCol[tid]=col;
	    mTowerId[ns][row][col]=tid;
	    tid++;
	  }else{
	    mTowerId[ns][row][col]=-1;
	  }
	}
      }
      for(int row=0; row<maxRowBlock(); row++){
	for(int col=0; col<maxColBlock(); col++){
	  if(doesBlockExist(ns,row,col)){
	    mBlockNS [bid]=ns;
	    mBlockRow[bid]=row;
	    mBlockCol[bid]=col;
	    mBlockId[ns][row][col]=bid;
	    //printf("%4d %2d %2d\n",bid,row,col);
	    bid++;
	  }else{
	    mBlockId[ns][row][col]=-1;
	  }
	}
      }    
      for(int row=0; row<maxRowFeebd(); row++){
	for(int col=0; col<maxColFeebd(); col++){
	  if(doesFeebdExist(ns,row,col)){
	    mFeebdNS[fid]=ns;
	    mFeebdRow[fid]=row;
	    mFeebdCol[fid]=col;
	    mFeebdId[ns][row][col]=fid;
	    fid++;
	  }else{
	    mFeebdId[ns][row][col]=-1;
	  }
	}
      }
    }
  }
};
