//=================================================================================
//C implementation of Fuzzy SLIC (Precise superpixel number control version) with Matlab Interface
//
//
/* Copyright (c) 2018, Chong WU All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
  
* Neither the name of City University of Hong Kong nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*Usage:
 * You need Matlab and a C/C++ compiler (latest version) to compile this code before using it.
 *
 * Input variables are:
 * 1. 8 bit images (color or grayscale)
 * 2. Number of intended superpixels (optional, default is 200)
 * 3. Compactness coefficient (optional, default is 10), used to control the color sensitive.
 * 4. Optional input (Use 1 to select Fuzzy SLICNC, use 0 to select Fuzzy SLIC, default is 0)
 * 5. Optional input (Amplification coefficient, default is 0.2)
 * Ouputs are:
 * 1. Labels (in raster scan order)
 * 2. Number of labels (final output superpixels) in the image
 *
 * Syntax:
 * [labels, numlabels] = FuzzySLIC(img,200,15);
 * [labels, numlabels] = FuzzySLIC(img,200,15,1);
 * [labels, numlabels] = FuzzySLIC(img,200,15,1,0.2);
 * NOTES:
 * Number of returned superpixels may be different from the input number of superpixels, if you select Fuzzy SLIC.
 * 
 * If you use this code please cite:

   C. Wu et al., "Fuzzy SLIC: Fuzzy Simple Linear Iterative Clustering," in IEEE Transactions on Circuits and Systems for Video Technology, doi: 10.1109/TCSVT.2020.3019109.

 */


#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>

double RGB2rgb(double val){
    if(val <= 0.04045){
        return val/12.92;
    }else{
        return pow((val+0.055)/1.055,2.4);
    }
}

void Rgb2Lab(int* rvec,int* gvec,int* bvec,int length,double* luminosityvec,
             double* channel_Avec, double* channel_Bvec){
    int index;
    int tempR,tempG,tempB;
    const double eps = 0.008856;
    const double kap = 903.3;
    const double XRef = 0.950456;
    const double YRef = 1.0;
    const double ZRef = 1.088754;
    const double coeX[3] = {0.4124564,0.3575761,0.1804375};
    const double coeY[3] = {0.2126729,0.7151522,0.0721750};
    const double coeZ[3] = {0.0193339,0.1191920,0.9503041};
    double X,Y,Z;
    double Red,Green,Blue;
    double red,green,blue;
    double xRef,yRef,zRef;
    double fx,fy,fz;
    double lValue,aValue,bValue;
    
    for (index = 0; index < length; ++index) {
        tempR = rvec[index];
        tempG = gvec[index];
        tempB = bvec[index];
        Red = tempR/255.0;
        Green = tempG/255.0;
        Blue = tempB/255.0;
        
        red = RGB2rgb(Red);
        green = RGB2rgb(Green);
        blue = RGB2rgb(Blue);
        
        X = red*coeX[0] + green*coeX[1] + blue*coeX[2];
        Y = red*coeY[0] + green*coeY[1] + blue*coeY[2];
        Z = red*coeZ[0] + green*coeZ[1] + blue*coeZ[2];
        
        xRef = X/XRef;
        yRef = Y/YRef;
        zRef = Z/ZRef;
        
        if(xRef>eps){fx = pow(xRef,1.0/3.0);}
        else{fx = (kap*xRef+16.0)/116.0;}
        if(yRef>eps){fy = pow(yRef,1.0/3.0);}
        else{fy = (kap*yRef+16.0)/116.0;}
        if(zRef>eps){fz = pow(zRef,1.0/3.0);}
        else{fz = (kap*zRef+16.0)/116.0;}
        
        lValue = 116.0*fy-16.0;
        aValue = 500.0*(fx-fy);
        bValue = 200.0*(fy-fz);
        
        luminosityvec[index] = lValue;
        channel_Avec[index] = aValue;
        channel_Bvec[index] = bValue;
        
        
    }
    
}

void getInitialClusterCentroids(int GridStep,int width,int height,
                                int* clusterlabels,int* numclusters){
    int temp,index;
    int xsteps,ysteps;
    int error1,error2;
    double error1perstep,error2perstep;
    int xthreshold,ythreshold;
    int axisX,axisY;
    int xerror1,yerror2;
    int clusterx,clustery;
    
    xsteps = (0.5+(double)(width)/(double)(GridStep));
    ysteps = (0.5+(double)(height)/(double)(GridStep));
    
    error1 = width - GridStep*xsteps;
    if(error1<0){
        xsteps--;
        error1 = width - GridStep*xsteps;
    }
    error2 = height - GridStep*ysteps;
    if(error2<0){
        ysteps--;
        error2 = height - GridStep*ysteps;
    }
    error1perstep = (double)(error1)/(double)(xsteps);
    error2perstep = (double)(error2)/(double)(ysteps);
    
    xthreshold = GridStep/2;
    ythreshold = GridStep/2;
    
    index = 0;
    for (axisY = 0;axisY<ysteps;++axisY){
        yerror2 = axisY*error2perstep;
        for(axisX = 0;axisX<xsteps;++axisX){
            xerror1 = axisX*error1perstep;
            clusterx = (axisX*GridStep+xthreshold+xerror1);
            clustery = (axisY*GridStep+ythreshold+yerror2);
            temp = clustery*width+clusterx;
            clusterlabels[index] = temp;
            ++index;
        }
    }
    *numclusters = index;
    
}

void FuzzySLIC(double* luminosityvec, double* channel_Avec, double* channel_Bvec, double* clusterl, double* clustera, double* clusterb, double* clusterx, double* clustery, int width, int height, int numclusters, int* oldlabels, int GridStep, double comp_coef)
{
    int x1, y1, x2, y2;
    double l, a, b, mid1,mid2,mid3,mid4,mid5,mid6;
    double dist;
    double distxy;
    int itr;
    int n;
    int x,y;
    int i,i2;
    int k,pos;
    double possi;
    int sz = width*height;
    const int numk = numclusters;
    int offset = GridStep;
    
    int GFSrange = 8;
    
    double* Hinv           = mxMalloc(sizeof(double)*numk);
    double* U1             = mxMalloc(sizeof(double)*numk);
    double* Ux             = mxMalloc(sizeof(double)*numk);
    double* Uy             = mxMalloc(sizeof(double)*numk);
    double* Ul             = mxMalloc(sizeof(double)*numk);
    double* Ua             = mxMalloc(sizeof(double)*numk);
    double* Ub             = mxMalloc(sizeof(double)*numk);
    
    double (*U)[GFSrange]         = mxMalloc(sizeof(double)*sz*GFSrange);
    double (*H)[GFSrange]         = mxMalloc(sizeof(double)*sz*GFSrange);
    int (*memb)[GFSrange]         = mxMalloc(sizeof(int)*sz*GFSrange);
    int* count             = mxMalloc(sizeof(int)*sz);
    double (*Lab_mat_c)[GFSrange] = mxMalloc(sizeof(double)*sz*GFSrange);
    int* x11 = mxMalloc(sizeof(int)*sz);
    int* y11 = mxMalloc(sizeof(int)*sz);
    
    double invwt = 1.0/((GridStep/comp_coef)*(GridStep/comp_coef));
    
    
    
    for( itr = 0; itr < 10; itr++ )
    {
        for(k = 0; k < sz; k++)
        {
            count[k] = 0;
            for (i = 0; i < GFSrange; i++)
            {
                /* code */
                Lab_mat_c[k][i] = 0;
                H[k][i] = 0;
                U[k][i] = 0;

            }
        }
        
        for(k = 0; k < numk; k++)
        {
            Hinv[k] = 0;
            U1[k] = 0;
            Ux[k] = 0;
            Uy[k] = 0;
            Ul[k] = 0;
            Ua[k] = 0;
            Ub[k] = 0;
            
        }
        for( n = 0; n < numk; n++ )
        {
            x1 = clusterx[n]-offset; if(x1 < 0) x1 = 0;
            y1 = clustery[n]-offset; if(y1 < 0) y1 = 0;
            x2 = clusterx[n]+offset; if(x2 > width)  x2 = width;
            y2 = clustery[n]+offset; if(y2 > height) y2 = height;
            
            for( y = y1; y < y2; y++ )
            {
                for( x = x1; x < x2; x++ )
                {
                    i = y*width + x;
                    x11[i] = x + 1;
                    y11[i] = y + 1;
                    l = luminosityvec[i];
                    a = channel_Avec[i];
                    b = channel_Bvec[i];
                    
                    dist =            (l - clusterl[n])*(l - clusterl[n]) +
                    (a - clustera[n])*(a - clustera[n]) +
                    (b - clusterb[n])*(b - clusterb[n]);
                    
                    distxy =        (x  - clusterx[n])*(x  - clusterx[n]) + (y  - clustery[n])*(y  - clustery[n]);
                    
                    dist += distxy*invwt;
                    
                    dist = sqrt(dist);
                    
                    //if(dist < 0.000000001)
                    //{
                        // if (count[i] < 3)
                        // {
                        //     if (count[i] == 0){
                        //         memb[i][count[i]] = n;
                        //         Lab_mat_c[i][count[i]] = 0.000000001;
                        //         count[i] += 1;
                        //     }
                        //     if (count[i] == 1){
                        //         memb[i][count[i]-1] = n;
                        //         Lab_mat_c[i][count[i]-1] = 0.000000001;
                        //     }
                        //     if (count[i] == 2){
                                
                        //         memb[i][count[i]-2] = n;
                        //         Lab_mat_c[i][count[i]-2] = 0.000000001;
                        //         count[i] -= 1;
                                
                        //     }
                            
                            
                        // }else{
                            
                        //     memb[i][0] = n;
                        //     Lab_mat_c[i][0] = 0.000000001;
                        //     memb[i][1] = n;
                        //     Lab_mat_c[i][1] = 0.000000001;
                        //     memb[i][2] = n;
                        //     Lab_mat_c[i][2] = 0.000000001;
                            
                        // }
                        
                        
                        
                    //}else{
                        if(count[i] < GFSrange)
                        {
                            memb[i][count[i]] = n;
                            Lab_mat_c[i][count[i]] = dist;
                            count[i] += 1;
                            
                        }else{


                            int MaxIndex = 0;
                            for(int index=MaxIndex+1;index<GFSrange;index++){
                                if (Lab_mat_c[i][index] > Lab_mat_c[i][MaxIndex]){
                                    MaxIndex = index;
                                }
                            }
                            pos = MaxIndex;
                            if (dist < Lab_mat_c[i][pos]) {
                                memb[i][pos] = n;
                                Lab_mat_c[i][pos] = dist;
                                
                            }else{
                                
                                
                            }
                            
                            
                            
                        }
                        
                    //}
                }
            }
        }
        
        for (i = 0; i < sz; i++){
            

                for (k = 0; k < count[i]; k++)
                {
                    mid1 = Lab_mat_c[i][k];
                    mid2 = 0;
                    for (int j = 0; j < count[i]; j++)
                    {
                        mid3 = 1.0/Lab_mat_c[i][j];
                        if (k == j){
                            mid2 = mid2 + 1.0;
                        }else
                        {
                            mid2 = mid2 + (mid1*mid3)*(mid1*mid3);
                        }
                        

                       
                    }
                    
                    U[i][k] = 1.0/mid2;
                    
                }

            
            
            
        }
        
        for (i = 0; i < sz; i++){
            
           
                
                if (((x11[i] - 1) >= 1) && ((y11[i] - 1) >= 1)){
                    i2 =  (y11[i] - 2)*width + x11[i] - 2;
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }                 
                }
                if (((x11[i] + 1) <= width) && ((y11[i] + 1) <= height)){
                    i2 =  (y11[i])*width + x11[i];
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }                   
                }
                if (((x11[i] + 1) <= width) && ((y11[i] - 1) >= 1)){
                    i2 =  (y11[i] - 2)*width + x11[i];
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }                     
                }
                if (((x11[i] - 1) >= 1) && ((y11[i] + 1) <= height)){
                    i2 =  (y11[i])*width + x11[i] - 2;
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }                     
                }
                if ((x11[i] - 1) >= 1){
                    i2 =  (y11[i] - 1)*width + x11[i] - 2;
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }  
                }
                if ((y11[i] - 1) >= 1){
                    i2 =  (y11[i] - 2)*width + x11[i] - 1;
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }  
                }
                if ((x11[i] + 1) <= width){
                    i2 =  (y11[i] - 1)*width + x11[i];
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                                
                            }
                            
                        }
                    }  
                }
                if ((y11[i] + 1) <= height){
                    i2 =  (y11[i])*width + x11[i] - 1;
                    for (int j = 0; j < count[i]; j++)
                    {
                        for (k = 0; k < count[i2]; k++)
                        {
                            if (memb[i2][k] == memb[i][j]){
                                    
                                H[i][j] = H[i][j] + U[i2][k];
                               
                            }
                            
                        }
                    }  
                }
                for (int j = 0; j < count[i]; j++)
                {
                    Hinv[memb[i][j]] = Hinv[memb[i][j]] + H[i][j]*H[i][j];
                }
            
            
        }
        
        for (i = 0; i < sz; i++){
            

            for (int j = 0; j < count[i]; j++)
            {
                U[i][j] = H[i][j]*H[i][j]/Hinv[memb[i][j]];
                U1[memb[i][j]] = U1[memb[i][j]] + U[i][j]*U[i][j];
                Ux[memb[i][j]] = Ux[memb[i][j]] + U[i][j]*U[i][j]*(x11[i]-1);
                Uy[memb[i][j]] = Uy[memb[i][j]] + U[i][j]*U[i][j]*(y11[i]-1);
                Ul[memb[i][j]] = Ul[memb[i][j]] + U[i][j]*U[i][j]*luminosityvec[i];
                Ua[memb[i][j]] = Ua[memb[i][j]] + U[i][j]*U[i][j]*channel_Avec[i];
                Ub[memb[i][j]] = Ub[memb[i][j]] + U[i][j]*U[i][j]*channel_Bvec[i];


            }
            int MaxIndex1 = 0;
            if (count[i]>1)
            {
                for(int index=MaxIndex1+1; index<count[i]; index++){
                    if (U[i][index] > U[i][MaxIndex1]){
                        MaxIndex1 = index;
                    }
                }          
            }
            
              
            oldlabels[i]  = memb[i][MaxIndex1];
            
            
            
        }
        
        
        
        for( k = 0; k < numk; k++ )
        {
            
            clusterl[k] = Ul[k]/U1[k];
            clustera[k] = Ua[k]/U1[k];
            clusterb[k] = Ub[k]/U1[k];
            clusterx[k] = Ux[k]/U1[k];
            clustery[k] = Uy[k]/U1[k];
            
        }
        
        
    }
    
    mxFree(Hinv);
    mxFree(H);
    mxFree(U);
    mxFree(U1);
    mxFree(Ul);
    mxFree(Ua);
    mxFree(Ub);
    mxFree(Ux);
    mxFree(Uy);
    mxFree(Lab_mat_c);
    mxFree(count);
    mxFree(y11);
    mxFree(x11);
    mxFree(memb);
    
}

void ConnectivityPostProcessing(int* preLabels,int width,int height,
                                int numSuperpixels,int* postLabels,int* numOfPostLabels){
    int i,j;
    int inner,outer,total;
    int x,y;
    int index;
    int oldIndex,newIndex,adjacentLabel;
    int numOfLabels;
    const int neighborx4[4] = {-1,0,1,0};
    const int neighbory4[4] = {0,-1,0,1};
    const int size = width*height;
    const int superpixelSize = size/numSuperpixels;
    int* xvec =  (int*)mxMalloc(sizeof(int)*superpixelSize*10);
    int* yvec =  (int*)mxMalloc(sizeof(int)*superpixelSize*10);
    
    for(i=0;i<size;++i){
        postLabels[i] = -1;
    }
    oldIndex = 0;
    adjacentLabel = 0;
    numOfLabels = 0;
    for(i = 0;i<height;++i){
        for(j = 0;j<width;++j){
            if(postLabels[oldIndex]<0){
                postLabels[oldIndex] = numOfLabels;
                xvec[0] = j;
                yvec[0] = i;
                
                for(inner = 0;inner<4;inner++){
                    x = xvec[0] + neighborx4[inner];
                    y = yvec[0] + neighbory4[inner];
                    if((x>=0&&x<width)&&(y>=0&&y<height)){
                        newIndex = y*width + x;
                        if(postLabels[newIndex]>=0){
                            adjacentLabel = postLabels[newIndex];
                        }
                    }
                }
                total = 1;
                for(outer = 0;outer<total;++outer){
                    for(inner = 0;inner<4;++inner){
                        x = xvec[outer] + neighborx4[inner];
                        y = yvec[outer] + neighbory4[inner];
                        if((x>=0&&x<width)&&(y>=0&&y<height)){
                            newIndex = y*width + x;
                            if(postLabels[newIndex]<0&&preLabels[oldIndex] == preLabels[newIndex]){
                                xvec[total] = x;
                                yvec[total] = y;
                                postLabels[newIndex] = numOfLabels;
                                total++;
                            }
                        }
                    }
                }
                if(total<=superpixelSize>>2){
                    for(outer = 0;outer<total;++outer){
                        index = yvec[outer]*width+xvec[outer];
                        postLabels[index] = adjacentLabel;
                    }
                    numOfLabels--;
                }
                numOfLabels++;
            }
            oldIndex++;
        }
    }
    *numOfPostLabels = numOfLabels;
    
    mxFree(xvec);
    mxFree(yvec);
}


void CombineSuperpixel(double* lvec, double* avec, double* bvec, int width, int height, int finalNumberOfLabels, int* clabels, int step, int exceptnum, double com)
{
    
    int i, j;
    int ind;
    int r, c;
    int k;
    int sz = width*height;
    int numk = finalNumberOfLabels;
    int diff = finalNumberOfLabels - exceptnum;
    int record = 1;
    double* clustersize = mxMalloc(sizeof(double)*numk);
    double* inv         = mxMalloc(sizeof(double)*numk);
    double* sigmal      = mxMalloc(sizeof(double)*numk);
    double* sigmaa      = mxMalloc(sizeof(double)*numk);
    double* sigmab      = mxMalloc(sizeof(double)*numk);
    double* sigmax      = mxMalloc(sizeof(double)*numk);
    double* sigmay      = mxMalloc(sizeof(double)*numk);
    double invwt = 1.0/((step/com)*(step/com));
    
    double* labelsize = mxMalloc(sizeof(double)*numk);
    double* labelindex = mxMalloc(sizeof(double)*numk);
//    printf("%d",diff);
//    printf(" ");
    //-----------------------------------------------------------------
    // Recalculate the centroid and store in the seed values
    //-----------------------------------------------------------------
    for(k = 0; k < numk; k++)
    {
        sigmal[k] = 0;
        sigmaa[k] = 0;
        sigmab[k] = 0;
        sigmax[k] = 0;
        sigmay[k] = 0;
        clustersize[k] = 0;
    }
    
    ind = 0;
    for( r = 0; r < height; r++ )
    {
        for( c = 0; c < width; c++ )
        {
            ind = r*width + c;

                sigmal[clabels[ind]] += lvec[ind];
                sigmaa[clabels[ind]] += avec[ind];
                sigmab[clabels[ind]] += bvec[ind];
                sigmax[clabels[ind]] += c;
                sigmay[clabels[ind]] += r;
                clustersize[clabels[ind]] += 1.0;
            
            
        }
    }
    
    for( k = 0; k < numk; k++ )
    {
        if( clustersize[k] <= 0 ) clustersize[k] = 1;
        inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
    }
    
    
    for( k = 0; k < numk; k++ )
    {
        sigmal[k] *= inv[k];
        sigmaa[k] *= inv[k];
        sigmab[k] *= inv[k];
        sigmax[k] *= inv[k];
        sigmay[k] *= inv[k];
    }
    
    for (k = 0; k < numk; k++)
    {
        labelindex[k] = k;
        labelsize[k] = clustersize[k];
    }
    
    for (i = 0; i < numk; i++)
    {
        for (j = 0; j < numk- i - 1; j++)
        {
            if (clustersize[j] > clustersize[j + 1])
            {
                int temp = clustersize[j];
                clustersize[j] = clustersize[j + 1];
                clustersize[j + 1] = temp;
                
                int ind_temp = labelindex[j];
                labelindex[j] = labelindex[j + 1];
                labelindex[j + 1] = ind_temp;
            }
        }
    }
    
    int labeldelate[diff];
    int nearindex[diff];
    for(i = 0; i < diff; i++)
    {
        labeldelate[i] = labelindex[i];
        nearindex[i] = 0;
    }
    
    int labelnear[diff][10];
    
    
    for(i = 0; i < diff; i++)
    {
        for(j = 0; j < 10; j++)
            labelnear[i][j] = labeldelate[i];
    }
    
    int now = 0, judge = 0, near1 = 0, near2 = 0, near3 = 0, near4 = 0, ind1 = 0, ind2 = 0, ind3 = 0, ind4 = 0, temps = 0;
    ind = 0;
    
    for( r = 0; r < height; r++ )
    {
        for( c = 0; c < width; c++ )
        {
            
            ind = r*width + c;
            judge = 0;
            now = clabels[ind];
            for(k = 0; k < diff; k++)
            {
                if(now == labeldelate[k])
                {
                    judge = 1;
                    break;
                }
            }
            
            if(judge == 1)
            {
                ind1 = ind - 1; if(ind1 < 0) ind1 = 0;
                ind2 = ind + 1; if(ind2 > (sz-1)) ind2 = sz-1;
                ind3 = ind - width; if(ind3 < 0) ind3 = ind;
                ind4 = ind + width; if(ind4 > (sz-1)) ind4 = ind;
                
                near1 = clabels[ind1];
                near2 = clabels[ind2];
                near3 = clabels[ind3];
                near4 = clabels[ind4];
                
                temps = 1;
                for(i = 0; i < diff; i++)
                {
                    if(near1 == labeldelate[i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near1;nearindex[k] += 1;}
                
                
                
                temps = 1;
                for(i = 0; i < diff; i++)
                {
                    if(near2 == labeldelate[i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near2;nearindex[k] += 1;}
                
                temps = 1;
                for(i = 0; i < diff; i++)
                {
                    if(near3 == labeldelate[i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near3;nearindex[k] += 1;}
                
                temps = 1;
                for(i = 0; i < diff; i++)
                {
                    if(near4 == labeldelate[i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near4;nearindex[k] += 1;}
                
                
            }
            
            
            
        }
    }
    
    for( r = 0; r < height; r++ )
    {
        for( c = 0; c < width; c++ )
        {
            
            ind = r*width + c;
            judge = 0;
            now = clabels[ind];
            for(k = 0; k < diff; k++)
            {
                if(now == labeldelate[k])
                {
                    judge = 1;
                    break;
                }
            }
            
            if(judge == 1)
            {
                ind1 = ind - 1; if(ind1 < 0) ind1 = 0;
                ind2 = ind + 1; if(ind2 > (sz-1)) ind2 = sz-1;
                ind3 = ind - width; if(ind3 < 0) ind3 = ind;
                ind4 = ind + width; if(ind4 > (sz-1)) ind4 = ind;
               
                near1 = clabels[ind1];
                near2 = clabels[ind2];
                near3 = clabels[ind3];
                near4 = clabels[ind4];
                
                temps = 1;
                for(i = 0; i < 10; i++)
                {
                    if(near1 == labelnear[k][i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near1;nearindex[k] += 1;}
                
                
                
                temps = 1;
                for(i = 0; i < 10; i++)
                {
                    if(near2 == labelnear[k][i])
                    {
                        temps = 0;
                        break;
                    }
                    
                }
                
                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near2;nearindex[k] += 1;}
                
                temps = 1;
                for(i = 0; i < 10; i++)
                {
                    if(near3 == labelnear[k][i])
                    {
                        temps = 0;
                        break;
                    }

                }

                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near3;nearindex[k] += 1;}

                temps = 1;
                for(i = 0; i < 10; i++)
                {
                    if(near4 == labelnear[k][i])
                    {
                        temps = 0;
                        break;
                    }

                }

                if(temps == 1 && nearindex[k]<10) {labelnear[k][nearindex[k]] = near4;nearindex[k] += 1;}

                
            }

            

        }
    }
    

    
    int labeltrans[diff];
    int labelnews[diff];
    for (i=0; i<diff; i++) {
        labeltrans[i] = -1;
        labelnews[i] = -1;
    }
    double distmin[diff];
    int dist, distxy, index, jud;
    near1 = 0;
    near2 = 0;
    
    
    
    while (record) {

        for (i = 0; i < diff; i++)
        {
            distmin[i] = DBL_MAX;
            index = labeldelate[i];
            ind = 0;
            if (labeltrans[i] < 0) {
                jud = 0;int tempjud;
                for (j = 0; j < 10; j++)
                {
                    tempjud = 1;
                    temps = labelnear[i][j];
                    near1 = temps;
                    for(k = 0; k < diff; k++)
                    {
                        if (temps == index) {
                            tempjud = 0;
                            jud++;
                            break;
                        }
                        
                        if (temps == labeldelate[k] && labelnews[k] < 0) {
                            tempjud = 0;
                            jud++;
                            break;
                        }
                        

                        if (temps == labeldelate[k] && labelnews[k] >= 0 && labeldelate[k] != index) {
                            temps = labelnews[k];
                            break;
                        }
                    }
                    
                    if (tempjud == 0) {
                        continue;
                    }
                    
         
                    
                    if (temps >= numk) continue;
                    if (temps < 0) continue;
                    dist =          (sigmal[index] - sigmal[temps])*(sigmal[index] - sigmal[temps]) +
                            (sigmaa[index] - sigmaa[temps])*(sigmaa[index] - sigmaa[temps]) +
                            (sigmab[index] - sigmab[temps])*(sigmab[index] - sigmab[temps]);
                    distxy =        (sigmax[index] - sigmax[temps])*(sigmax[index] - sigmax[temps]) +
                            (sigmay[index] - sigmay[temps])*(sigmay[index] - sigmay[temps]);
                    
                    dist += distxy*invwt;
                    if(dist < distmin[i])
                    {
                        distmin[i] = dist;
                        labeltrans[i] = temps;
             
                        labelnews[i] = temps;
                    }
                }
            }
        }
        
        for(i = 0; i<diff; i++){
            record = 0;
            if (labeltrans[i] < 0) {
                record = 1;
                break;
            }
            
        }
        
    }
        
        ind = 0;
        for( r = 0; r < height; r++ )
        {
            for( c = 0; c < width; c++ )
            {
                judge = 0;
                ind = r*width + c;
                now = clabels[ind];
                
                for(k = 0; k < diff; k++)
                {
                    if(now == labeldelate[k])
                    {
                        judge = 1;
                        break;
                    }
                }
                
                if(judge == 1 && labeltrans[k] >= 0)
                {
                    clabels[ind] = labeltrans[k];
                }
            }
        }
        
    
    
//    for (i = 0; i<diff; i++) {
//        printf("%d",labeltrans[i]);
//        printf(" ");
//    }
//    for (i = 0; i<diff; i++) {
//        printf("%d",labeldelate[i]);
//        printf(" ");
//    }
//    for(i = 0; i < diff; i++)
//    {
//        for(j = 0; j < 10; j++)
//        {
//        printf("%d",labelnear[i][j]);
//        printf(" ");
//        }
//    }
    
    mxFree(sigmal);
    mxFree(sigmaa);
    mxFree(sigmab);
    mxFree(sigmax);
    mxFree(sigmay);
    mxFree(clustersize);
    mxFree(inv);
    mxFree(labelsize);
    mxFree(labelindex);
    
    
}


void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    int width;
    int height;
    int size;
    int indexOfVec,indexOfImg,i;
    int x,y;
    int *rvec,*gvec,*bvec;
    int* oldlabels;
    int* newlabels;
    double *luminosityvec,*channel_Avec,*channel_Bvec;
    int GridStep;
    int* clusterlabels;
    int numclusters;
    double *clusterx,*clustery;
    double *clusterl,*clustera,*clusterb;
    const mwSize* dimensions;
    int* outputNumOfSuperpixels;
    int* ouputLabels;
    int numOfPostLabels;
    unsigned char* imgbytes;
    int numParameters;
    int numOfSuperpixels = 200;
    double comp_coef = 10;
    
    int numcontrol = 0;
    double num_coef = 0.2;
    
    
    int exceptnum, exceptnum1;
    
    if(nrhs<1){
        mexErrMsgTxt("At least one argument is required.") ;
    } else if (nrhs>5){
        mexErrMsgTxt("Too many input arguments.");
    }
    if(nlhs!=2){
        mexErrMsgIdAndTxt("Fuzzy SLIC:nlhs","Two outputs required, a labels and the number of labels, i.e superpixels.");
    }
    
    numParameters = mxGetNumberOfElements(prhs[0]);
    mwSize numOfDims = mxGetNumberOfDimensions(prhs[0]);
    dimensions = mxGetDimensions(prhs[0]);
    imgbytes = (unsigned char*)mxGetData(prhs[0]);
    height = dimensions[0];
    width = dimensions[1];
    size = width*height;
    
    numOfSuperpixels = mxGetScalar(prhs[1]);
    comp_coef = mxGetScalar(prhs[2]);
    if (nrhs==4) {
        numcontrol = mxGetScalar(prhs[3]);
        
    }
    if (nrhs==5) {
        numcontrol = mxGetScalar(prhs[3]);
        num_coef = mxGetScalar(prhs[4]);
    }
    
    exceptnum1 = numOfSuperpixels;
    
    
    rvec = (int*)mxMalloc(sizeof(int)*size);
    gvec = (int*)mxMalloc(sizeof(int)*size);
    bvec = (int*)mxMalloc(sizeof(int)*size);
    
    luminosityvec = (double*)mxMalloc(sizeof(double)*size);
    channel_Avec = (double*)mxMalloc(sizeof(double)*size);
    channel_Bvec = (double*)mxMalloc(sizeof(double)*size);
    oldlabels = (int*)mxMalloc(sizeof(int)*size);
    newlabels = (int*)mxMalloc(sizeof(int)*size);
    clusterlabels = (int*)mxMalloc(sizeof(int)*size);
    
    if(numParameters/size == 1){
        for(x = 0,indexOfImg = 0;x<width;++x){
            for(y = 0;y<height;++y){
                indexOfVec = y*width+x;
                luminosityvec[indexOfVec] = imgbytes[indexOfImg];
                channel_Avec[indexOfVec] = imgbytes[indexOfImg];
                channel_Bvec[indexOfVec] = imgbytes[indexOfImg];
                indexOfImg++;
            }
        }
    }
    else{
        for(x = 0,indexOfImg = 0;x<width;++x){
            for(y = 0;y<height;++y){
                indexOfVec = y*width+x;
                rvec[indexOfVec] = imgbytes[indexOfImg];
                gvec[indexOfVec] = imgbytes[indexOfImg+size];
                bvec[indexOfVec] = imgbytes[indexOfImg+size+size];
                ++indexOfImg;
            }
        }
        Rgb2Lab(rvec,gvec,bvec,size,luminosityvec,channel_Avec,channel_Bvec);
    }
    
    GridStep = sqrt((double)(size)/(double)(numOfSuperpixels))+0.5;
    getInitialClusterCentroids(GridStep,width,height,clusterlabels,&numclusters);
    
    exceptnum = numclusters;
    
    
    clusterx = (double*)mxMalloc(sizeof(double)*numclusters);
    clustery = (double*)mxMalloc(sizeof(double)*numclusters);
    clusterl = (double*)mxMalloc(sizeof(double)*numclusters);
    clustera = (double*)mxMalloc(sizeof(double)*numclusters);
    clusterb = (double*)mxMalloc(sizeof(double)*numclusters);
    
    for(i = 0;i<numclusters;++i){
        clusterx[i] = clusterlabels[i]%width;
        clustery[i] = clusterlabels[i]/width;
        clusterl[i] = luminosityvec[clusterlabels[i]];
        clustera[i] = channel_Avec[clusterlabels[i]];
        clusterb[i] = channel_Bvec[clusterlabels[i]];
    }
    
    FuzzySLIC(luminosityvec,channel_Avec,channel_Bvec,clusterl,clustera,clusterb,clusterx,
              clustery,width,height,numclusters,oldlabels,GridStep,comp_coef);
    ConnectivityPostProcessing(oldlabels,width,height,numOfSuperpixels,
                               newlabels,&numOfPostLabels);
     //printf("%d",numOfPostLabels);
    // Compute the loss-rate of superpixels
    
    mxFree(clusterx);
    mxFree(clustery);
    mxFree(clusterl);
    mxFree(clustera);
    mxFree(clusterb);
    
    if (numOfPostLabels < exceptnum && numcontrol == 1)
    {
        numOfSuperpixels = (numOfSuperpixels * exceptnum1) / numOfPostLabels + exceptnum1 * num_coef;
        
        
        // Recalculate
        GridStep = sqrt((double)(size)/(double)(numOfSuperpixels))+0.5;
        getInitialClusterCentroids(GridStep,width,height,clusterlabels,&numclusters);
        //printf("%d",numclusters);
        
        clusterx = (double*)mxMalloc(sizeof(double)*numclusters);
        clustery = (double*)mxMalloc(sizeof(double)*numclusters);
        clusterl = (double*)mxMalloc(sizeof(double)*numclusters);
        clustera = (double*)mxMalloc(sizeof(double)*numclusters);
        clusterb = (double*)mxMalloc(sizeof(double)*numclusters);
        
        for(i = 0;i<numclusters;++i){
            clusterx[i] = clusterlabels[i]%width;
            clustery[i] = clusterlabels[i]/width;
            clusterl[i] = luminosityvec[clusterlabels[i]];
            clustera[i] = channel_Avec[clusterlabels[i]];
            clusterb[i] = channel_Bvec[clusterlabels[i]];
        }
        
        FuzzySLIC(luminosityvec,channel_Avec,channel_Bvec,clusterl,clustera,clusterb,clusterx,
                  clustery,width,height,numclusters,oldlabels,GridStep,comp_coef);
        ConnectivityPostProcessing(oldlabels,width,height,numOfSuperpixels,
                                   newlabels,&numOfPostLabels);
        
        
        if (numOfPostLabels > exceptnum)
        {
            CombineSuperpixel(luminosityvec, channel_Avec, channel_Bvec, width, height, numOfPostLabels, newlabels, GridStep, exceptnum, comp_coef);
            numOfPostLabels = exceptnum;
        }
        mxFree(clusterx);
        mxFree(clustery);
        mxFree(clusterl);
        mxFree(clustera);
        mxFree(clusterb);
    }
    
    else
    {
        if (numOfPostLabels > exceptnum)
        {
            CombineSuperpixel(luminosityvec, channel_Avec, channel_Bvec, width, height, numOfPostLabels, newlabels, GridStep, exceptnum, comp_coef);
            numOfPostLabels = exceptnum;
        }
        
    }
    
    plhs[0] = mxCreateNumericMatrix(height,width,mxINT32_CLASS,mxREAL);
    ouputLabels = (int*)mxGetData(plhs[0]);
    for (x = 0,indexOfImg = 0;x<width;++x){
        for(y = 0;y<height;++y){
            indexOfVec = y*width+x;
            ouputLabels[indexOfImg] = newlabels[indexOfVec];
            indexOfImg++;
        }
    }
    
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    outputNumOfSuperpixels = (int*)mxGetData(plhs[1]);
    *outputNumOfSuperpixels = numOfPostLabels;
    
    
    
    
    mxFree(rvec);
    mxFree(gvec);
    mxFree(bvec);
    mxFree(luminosityvec);
    mxFree(channel_Avec);
    mxFree(channel_Bvec);
    mxFree(oldlabels);
    mxFree(newlabels);
    mxFree(clusterlabels);

}
