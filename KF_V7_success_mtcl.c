/*
 * @Author: your name
 * @Date: 2020-05-29 14:05:55
 * @LastEditTime: 2020-05-29 22:33:04
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: \testKF_V1\KF_V5.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TarDims 2
#define MeasDims 2
#define VectorDims 2 * TarDims
#define simTime 100
#define MtclTime 100
#define T 1.0
#define PI (double)(3.141592654)

/********************************类型定义及函数原型***************************************************/
typedef struct
{
    int row, col;
    double **mat;
} NewMat;

typedef struct RecordMatrix
{
    int row, col, page;
    double ***mat;
} NewMatRec;

void IniNewMat(NewMat *ResMat, int rowNew, int colNew);
void FreeMat(NewMat *ResMat);
void PrintMat(NewMat *PrintMatrix);
void SetF(NewMat *F);
void SetQ(NewMat *Q, double q);
void SetH(NewMat *H);
void SetR(NewMat *R, double *StandDev);
void MatSum(NewMat *tmp, NewMat *A, NewMat *B, int JudgeNum);
void MatMulti(NewMat *tmp, NewMat *A, NewMat *B);
// void MatInverse(NewMat *tmp, NewMat *A);
void MatTranspose(NewMat *tmp, NewMat *A);
double getA(NewMat *Mat, int n);
void getAStart(NewMat *Mat, int n, NewMat *AStar);
void getAInverse(NewMat *Res, NewMat *Mat);
void SetUniMat(NewMat *UnitMat);
void MatCopy(NewMat *CopyMat, NewMat *OriMat);
double gaussianRand(double mean, double variance);
void IniNewMatRec(NewMatRec *IniMat, int rowNew, int colNew, int pageNew);
void FreeMatRec(NewMatRec *ResMat);

/********************************主**程**序***************************************************/
int main(int argc, char const *argv[])
{

    clock_t start, finish;
    double duration; /* 测量一dao个事件持续专的属时间*/

    // printf("%d", RAND_MAX);
    //设置状态转移矩阵
    NewMat F;
    IniNewMat(&F, 2 * TarDims, 2 * TarDims);
    SetF(&F);
    NewMat F_tran;
    IniNewMat(&F_tran, F.col, F.row);
    MatTranspose(&F_tran, &F);
    // FreeMat(&F);
    // PrintMat(&F);

    // 设置过程噪声误差协方差矩阵
    NewMat Q;
    IniNewMat(&Q, 2 * TarDims, 2 * TarDims);
    double q = 0.01;
    SetQ(&Q, q);
    // PrintMat(&Q);
    // printf("%f",Q.mat[1][1]);

    // 设置传感器观测矩阵
    NewMat H;
    IniNewMat(&H, MeasDims, 2 * TarDims);
    SetH(&H);
    NewMat H_tran;
    IniNewMat(&H_tran, H.col, H.row);
    MatTranspose(&H_tran, &H);
    // PrintMat(&H);

    //设置传感器观测误差协方差矩阵
    NewMat R;
    double StandDev[MeasDims] = {10, 10}; //设置观测误差一倍标准差，单位：m
    IniNewMat(&R, MeasDims, MeasDims);
    SetR(&R, &StandDev[0]);
    // PrintMat(&R);

    // 设置状态估计均值矩阵--预测
    NewMat Xp;
    IniNewMat(&Xp, 2 * TarDims, 1);

    // 设置测量预测
    NewMat Zp;
    IniNewMat(&Zp, MeasDims, 1);

    // 设置状态估计误差协方差矩阵--预测
    NewMat Pp;
    IniNewMat(&Pp, 2 * TarDims, 2 * TarDims);

    // 设置状态估计均值矩阵--更新
    NewMat Xe;
    IniNewMat(&Xe, 2 * TarDims, 1);
    // PrintMat(&Xe);

    // 设置状态估计误差协方差矩阵--更新
    NewMat Pe;
    IniNewMat(&Pe, 2 * TarDims, 2 * TarDims);

    // 设置新息协方差矩阵
    NewMat Sk;
    IniNewMat(&Sk, MeasDims, MeasDims);

    // 设置卡尔曼增益矩阵
    NewMat K;
    IniNewMat(&K, 2 * TarDims, MeasDims);

    // 设置测量矩阵
    NewMat Zk;
    IniNewMat(&Zk, MeasDims, 1);

    // 设置单位阵
    NewMat UnitMat;
    IniNewMat(&UnitMat, 2 * TarDims, 2 * TarDims);
    SetUniMat(&UnitMat);

    // 设置记录矩阵
    NewMatRec X_true_Rec;
    NewMatRec Z_true_Rec;
    NewMatRec Xe_true_Rec;
    NewMatRec State_Esti_Diff;
    NewMatRec Staet_Obse_Diff;
    IniNewMatRec(&X_true_Rec, 2 * TarDims, simTime, MtclTime);
    IniNewMatRec(&Z_true_Rec, MeasDims, simTime, MtclTime);
    IniNewMatRec(&Xe_true_Rec, 2 * TarDims, simTime, MtclTime);
    IniNewMatRec(&X_true_Rec, 2 * TarDims, simTime, MtclTime);
    IniNewMatRec(&Staet_Obse_Diff, MeasDims, simTime, MtclTime);
    IniNewMatRec(&State_Esti_Diff, 2 * TarDims, simTime, MtclTime);
    // 设置记录矩阵
    // NewMat X_true_Rec;
    // IniNewMat(&X_true_Rec, 2 * TarDims, simTime);
    // NewMat X_esti_Rec;
    // IniNewMat(&X_esti_Rec, 2 * TarDims, simTime);
    // NewMat DiffRec;
    // IniNewMat(&DiffRec, 2*TarDims, simTime);
    // 设置相关参数
    double Vmax = 20;
    double X0[2 * TarDims] = {10, 10, 10, 10};

    srand((unsigned)time(NULL)); // 设置随机数种子--->这样才会产生不同的随机数，否则使用相同种子产生的随机数都是完全相同的！！！
    // 正常时域循环迭代
    for (int mtclInd = 0; mtclInd < MtclTime; mtclInd++)
    {
        printf("mtclInd = %d \n", mtclInd);
        start = clock();
        // 各种初始化
        NewMat X_True;
        IniNewMat(&X_True, 2 * TarDims, 1);
        X_True.mat[0][0] = X0[0];
        X_True.mat[1][0] = X0[1];
        X_True.mat[2][0] = X0[2];
        X_True.mat[3][0] = X0[3];
        NewMat X_True_WithQ;
        IniNewMat(&X_True_WithQ, 2 * TarDims, 1);

        //生成第一帧测量
        NewMat Z_true;
        IniNewMat(&Z_true, MeasDims, 1);
        MatMulti(&Z_true, &H, &X_True);
        // 用第一帧测量初始化
        Xe.mat[0][0] = Z_true.mat[0][0];
        Xe.mat[1][0] = Z_true.mat[1][0];
        Xe.mat[2][0] = 0;
        Xe.mat[3][0] = 0;
        Pe.mat[0][0] = R.mat[0][0];
        Pe.mat[1][1] = R.mat[1][1];
        Pe.mat[2][2] = pow(Vmax, 2) / ((double)TarDims + 2);
        Pe.mat[3][3] = pow(Vmax, 2) / ((double)TarDims + 2);
        // PrintMat(&R);
        // // printf("\n");
        // printf("X_True:");
        // PrintMat(&X_True);
        // printf("Z_true:");
        // PrintMat(&Z_true);
        // printf("初始化状态估计均值  Xe:\n");
        // PrintMat(&Xe);
        // printf("初始化状态估计误差协方差  Pe:\n");
        // PrintMat(&Pe);
        // // printf("\n");
        for (int simScan = 1; simScan < simTime; simScan++)
        {
            // printf("simScan = %d\n", simScan);
            NewMat tmpMatrix1;
            NewMat tmpMatrix2;
            NewMat tmpMatrix3;

            // 生成测量数据
            IniNewMat(&tmpMatrix1, X_True.row, X_True.col);
            MatCopy(&tmpMatrix1, &X_True);
            MatMulti(&X_True, &F, &tmpMatrix1);
            double RandomQ;
            for (int i = 0; i < 2 * TarDims; i++)
            {
                RandomQ = gaussianRand(0.0, Q.mat[i][i]);
                X_True_WithQ.mat[i][0] = X_True.mat[i][0] + RandomQ;
            }
            // printf("Q = %f\n", gaussianRand(0.0, Q.mat[0][0]));

            MatMulti(&Zk, &H, &X_True_WithQ);
            double RandomR;
            for (int j = 0; j < MeasDims; j++)
            {
                RandomR = gaussianRand(0.0, R.mat[j][j]);
                Zk.mat[j][0] += RandomR;
            }

            // printf("R = %f\n", gaussianRand(0.0, R.mat[0][0]));
            FreeMat(&tmpMatrix1);

            // 状态预测
            IniNewMat(&tmpMatrix1, F.row, Pe.col);
            MatMulti(&Xp, &F, &Xe); //Xp
            MatMulti(&tmpMatrix1, &F, &Pe);
            IniNewMat(&tmpMatrix2, tmpMatrix1.row, F_tran.col);
            MatMulti(&tmpMatrix2, &tmpMatrix1, &F_tran);
            MatSum(&Pp, &tmpMatrix2, &Q, 1); //Pp
            FreeMat(&tmpMatrix1);
            FreeMat(&tmpMatrix2);

            // 新息协方差
            IniNewMat(&tmpMatrix1, H.row, Pp.col);
            MatMulti(&tmpMatrix1, &H, &Pp);
            IniNewMat(&tmpMatrix2, tmpMatrix1.row, H_tran.col);
            MatMulti(&tmpMatrix2, &tmpMatrix1, &H_tran);
            MatSum(&Sk, &tmpMatrix2, &R, 1);

            FreeMat(&tmpMatrix1);
            FreeMat(&tmpMatrix2);

            // 卡尔曼增益矩阵
            IniNewMat(&tmpMatrix1, Pp.row, H_tran.col);
            MatMulti(&tmpMatrix1, &Pp, &H_tran);
            IniNewMat(&tmpMatrix2, Sk.row, Sk.col);
            getAInverse(&tmpMatrix2, &Sk);
            MatMulti(&K, &tmpMatrix1, &tmpMatrix2);
            // printf("inv(Sk): \n");
            // PrintMat(&tmpMatrix2);
            IniNewMat(&tmpMatrix3, Sk.row, Sk.col);
            MatMulti(&tmpMatrix3, &Sk, &tmpMatrix2);

            FreeMat(&tmpMatrix1);
            FreeMat(&tmpMatrix2);

            // 状态更新
            // 更新均值
            IniNewMat(&tmpMatrix1, Zp.row, Zp.col);
            MatMulti(&Zp, &H, &Xp);
            MatSum(&tmpMatrix1, &Zk, &Zp, -1);

            IniNewMat(&tmpMatrix2, K.row, tmpMatrix1.col);
            MatMulti(&tmpMatrix2, &K, &tmpMatrix1);
            // printf("K*残差: \n");
            // PrintMat(&tmpMatrix2);
            MatSum(&Xe, &tmpMatrix2, &Xp, 1);

            FreeMat(&tmpMatrix1);
            FreeMat(&tmpMatrix2);

            // 更新估计误差协方差矩阵
            IniNewMat(&tmpMatrix1, K.row, H.col);
            IniNewMat(&tmpMatrix2, UnitMat.row, tmpMatrix1.col);
            MatMulti(&tmpMatrix1, &K, &H);
            MatSum(&tmpMatrix2, &UnitMat, &tmpMatrix1, -1);
            MatMulti(&Pe, &tmpMatrix2, &Pp);

            // 输出估计误差
            // for (int i = 0; i < 2 * TarDims; i++)
            // {
            //     printf("%.3f  ", X_True_WithQ.mat[i][0] - Xe.mat[i][0]);
            // }
            // printf("\n");

            // 根据跟踪结果统计状态轨迹差值
            for (int calInd = 0; calInd < 2 * TarDims; calInd++)
            {
                X_true_Rec.mat[calInd][simScan][mtclInd] = X_True_WithQ.mat[calInd][0];
                Xe_true_Rec.mat[calInd][simScan][mtclInd] = Xe.mat[calInd][0];
                State_Esti_Diff.mat[calInd][simScan][mtclInd] = Xe_true_Rec.mat[calInd][simScan][mtclInd] - X_true_Rec.mat[calInd][simScan][mtclInd];
                if (calInd < TarDims)
                {
                    Z_true_Rec.mat[calInd][simScan][mtclInd] = Zk.mat[calInd][0];
                    Staet_Obse_Diff.mat[calInd][simScan][mtclInd] = Z_true_Rec.mat[calInd][simScan][mtclInd] - X_true_Rec.mat[calInd][simScan][mtclInd];
                }
            }

            // printf("X_true_WithQ: \n");
            // PrintMat(&X_True_WithQ);
            // printf("Zk: \n");
            // PrintMat(&Zk);
            // printf("Xp: \n");
            // PrintMat(&Xp);
            // printf("Sk: \n");
            // PrintMat(&Sk);
            // printf("Sk * inv(Sk): \n");
            // PrintMat(&tmpMatrix3);
            // printf("K: \n");
            // PrintMat(&K);
            // printf("Zp：\n");
            // PrintMat(&Zp);
            // printf("残差：\n");
            // PrintMat(&tmpMatrix1);
            // printf("Xe:\n");
            // PrintMat(&Xe);
            // printf("X_true_WithQ: \n");
            // PrintMat(&X_True_WithQ);

            // FreeMat(&tmpMatrix1);
            // FreeMat(&tmpMatrix2);
            // FreeMat(&tmpMatrix3);

            // printf("目标真实状态X_true_WithQ:\n");
            // PrintMat(&X_True_WithQ);
            // printf("目标真实测量Zk:\n");
            // PrintMat(&Zk);

            // MatMulti(&Zp, &H, &Xp);
            // printf("状态预测Xp:\n");
            // PrintMat(&Xp);
            // printf("状态预测协方差Pp:\n");
            // PrintMat(&Pp);

            // printf("新息协方差Sk:\n");
            // PrintMat(&Sk);

            // printf("目标状态估计更新Xe:\n");
            // PrintMat(&Xe);
            // printf("目标状态估计协方差更新Pe:\n");
            // PrintMat(&Pe);
        }
        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        printf("%f seconds\n", duration);
    }

    // 所有蒙特卡洛实验结束以后统计性能指标
    NewMatRec State_Esti_Diff_Square;
    IniNewMatRec(&State_Esti_Diff_Square, State_Esti_Diff.row, State_Esti_Diff.col, State_Esti_Diff.page);
    NewMatRec State_Obse_Diff_Square;
    IniNewMatRec(&State_Obse_Diff_Square, Staet_Obse_Diff.row, Staet_Obse_Diff.col, Staet_Obse_Diff.page);
    for (int cal_row = 0; cal_row < State_Esti_Diff.row; cal_row++)
    {
        for (int cal_col = 0; cal_col < State_Esti_Diff.col; cal_col++)
        {
            for (int cal_page = 0; cal_page < State_Esti_Diff.page; cal_page++)
            {
                State_Esti_Diff_Square.mat[cal_row][cal_col][cal_page] = pow(State_Esti_Diff.mat[cal_row][cal_col][cal_page], 2);
                if (cal_row < TarDims)
                {
                    State_Obse_Diff_Square.mat[cal_row][cal_col][cal_page] = pow(Staet_Obse_Diff.mat[cal_row][cal_col][cal_page], 2);
                }
            }
        }
    }

    // for (int i = 0; i < simTime; i++)
    // {
    //     printf("PosSquare = %.3f ", State_Obse_Diff_Square.mat[0][i][0]);
    //     printf("MeaSquare = %.3f ", State_Esti_Diff_Square.mat[0][i][0]);
    //     printf("\n");
    // }

    NewMatRec State_Esti_Sum;
    NewMatRec State_Obse_Sum;
    IniNewMatRec(&State_Esti_Sum, TarDims, simTime, MtclTime);
    IniNewMatRec(&State_Obse_Sum, 1, simTime, MtclTime);

    for (int cal_col = 0; cal_col < State_Esti_Diff.col; cal_col++)
    {
        for (int cal_page = 0; cal_page < State_Esti_Diff.page; cal_page++)
        {
            State_Esti_Sum.mat[0][cal_col][cal_page] = State_Esti_Diff_Square.mat[0][cal_col][cal_page] + State_Esti_Diff_Square.mat[1][cal_col][cal_page];
            State_Esti_Sum.mat[1][cal_col][cal_page] = State_Esti_Diff_Square.mat[2][cal_col][cal_page] + State_Esti_Diff_Square.mat[3][cal_col][cal_page];
            State_Obse_Sum.mat[0][cal_col][cal_page] = State_Obse_Diff_Square.mat[0][cal_col][cal_page] + State_Obse_Diff_Square.mat[1][cal_col][cal_page];
        }
    }

    // for (int i = 0; i < simTime; i++)
    // {
    //     printf("PosSquare = %.3f ", State_Esti_Sum.mat[0][i][0]);
    //     printf("MeaSquare = %.3f ", State_Obse_Sum.mat[0][i][0]);
    //     printf("\n");
    // }

    for (int cal_col = 0; cal_col < simTime; cal_col++)
    {
        for (int cal_page = 1; cal_page < MtclTime; cal_page++)
        {
            State_Esti_Sum.mat[0][cal_col][0] += State_Esti_Sum.mat[0][cal_col][cal_page];
            State_Esti_Sum.mat[1][cal_col][0] += State_Esti_Sum.mat[1][cal_col][cal_page];
            State_Obse_Sum.mat[0][cal_col][0] += State_Obse_Sum.mat[0][cal_col][cal_page];
        }
        State_Esti_Sum.mat[0][cal_col][0] /= MtclTime;
        State_Esti_Sum.mat[1][cal_col][0] /= MtclTime;
        State_Obse_Sum.mat[0][cal_col][0] /= MtclTime;
    }

    // for (int i = 0; i < simTime; i++)
    // {
    //     printf("PosSquare = %.3f ", State_Esti_Sum.mat[0][i][0]);
    //     printf("MeaSquare = %.3f ", State_Obse_Sum.mat[0][i][0]);
    //     printf("\n");
    // }

    NewMat RMSE_Pos, RMSE_Vel, RMSE_Pos_Z;
    IniNewMat(&RMSE_Pos, 1, simTime);
    IniNewMat(&RMSE_Vel, 1, simTime);
    IniNewMat(&RMSE_Pos_Z, 1, simTime);

    for (int cal_col = 0; cal_col < simTime; cal_col++)
    {
        RMSE_Pos.mat[0][cal_col] = sqrt(State_Esti_Sum.mat[0][cal_col][0]);
        RMSE_Vel.mat[0][cal_col] = sqrt(State_Esti_Sum.mat[1][cal_col][0]);
        RMSE_Pos_Z.mat[0][cal_col] = sqrt(State_Obse_Sum.mat[0][cal_col][0]);
        // printf("RMSE_Pos = %.3f    RMSE_Pos_Z = %.3f", RMSE_Pos.mat[0][cal_col], RMSE_Pos_Z.mat[0][cal_col]);
        if ((RMSE_Pos.mat[0][cal_col] - RMSE_Pos_Z.mat[0][cal_col]) > 0)
        {
            printf("无效 \n");
        }else
        {
            printf("有效 \n");
        }
    }

    // 释放相应内存
    FreeMat(&F);
    FreeMat(&Q);
    FreeMat(&H);
    FreeMat(&R);
    FreeMat(&Xp);
    FreeMat(&Pp);
    FreeMat(&Xe);
    FreeMat(&Pe);
    FreeMatRec(&X_true_Rec);
    FreeMatRec(&Xe_true_Rec);
    FreeMatRec(&Z_true_Rec);
    FreeMatRec(&Staet_Obse_Diff);
    FreeMatRec(&State_Esti_Diff);
    FreeMatRec(&State_Obse_Diff_Square);
    FreeMatRec(&State_Esti_Diff_Square);
    FreeMatRec(&State_Obse_Sum);
    FreeMatRec(&State_Esti_Sum);
    FreeMat(&RMSE_Pos);
    FreeMat(&RMSE_Vel);
    FreeMat(&RMSE_Pos_Z);
    // V1-报错--warning: attempt to free a non-heap object 'F''Q''H''R''Xp''Pp''Xe''Pe'
    return 0;
}

/********************************定**义**函**数***************************************************/
// 定义函数：初始化矩阵结构体
void IniNewMat(NewMat *ResMat, int rowNew, int colNew)
{
    // ResMat = (NewMat *)malloc(sizeof(NewMat));
    ResMat->row = rowNew;
    ResMat->col = colNew;
    ResMat->mat = (double **)malloc(rowNew * sizeof(double *));
    for (int i = 0; i < rowNew; i++)
    {
        ResMat->mat[i] = (double *)malloc(colNew * sizeof(double));
        if (ResMat->mat[i] == NULL)
        {
            puts("Memory allocations failed.Goodbye.");
            exit(EXIT_FAILURE);
        }
        else
        {
            // puts("Memory allocations success.Yeah.");
            for (int j = 0; j < colNew; j++)
                ResMat->mat[i][j] = 0.0;
        }
    }
    // PrintMat(ResMat);
}

// 定义函数：打印矩阵
void PrintMat(NewMat *PrintMatrix)
{
    for (int i = 0; i < PrintMatrix->row; i++)
    {
        for (int j = 0; j < PrintMatrix->col; j++)
            printf("%.10f  ", PrintMatrix->mat[i][j]);
        printf("\n");
    }
    printf("\n");
}

// 定义函数：设定状态转移矩阵
void SetF(NewMat *F)
{
    for (int i = 0; i < 2 * TarDims; i++)
    {
        F->mat[i][i] = 1;
        if (i < TarDims)
            F->mat[i][i + TarDims] = T;
    }
}

// 定义函数：设置过程噪声误差协方差矩阵
void SetQ(NewMat *Q, double q)
{
    for (int i = 0; i < TarDims; i++)
    {
        Q->mat[i][i] = q * pow(T, 3.0) / T;
        Q->mat[i][i + TarDims] = q * pow(T, 2.0) / T;
        Q->mat[i + TarDims][i] = q * pow(T, 2.0) / T;
        Q->mat[i + TarDims][i + TarDims] = q * T;
    }
}

// 定义函数：设置传感器观测矩阵-->线性测量
void SetH(NewMat *H)
{
    for (int i = 0; i < MeasDims; i++)
        H->mat[i][i] = 1;
}

// 定义函数：设置传感器观测噪声误差协方差矩阵
void SetR(NewMat *R, double *StandDev)
{
    for (int i = 0; i < MeasDims; i++)
        R->mat[i][i] = pow(StandDev[i], 2);
}

// 定义函数：设置单位矩阵
void SetUniMat(NewMat *UnitMat)
{
    for (int i = 0; i < 2 * TarDims; i++)
        UnitMat->mat[i][i] = 1;
}

// 定义函数：释放矩阵内存
void FreeMat(NewMat *ResMat)
{
    // free(ResMat->row);
    // free(ResMat->col);
    for (int i = 0; i < ResMat->row; i++)
        free(ResMat->mat[i]);
    // puts("free Memory Done.");
    free(ResMat->mat);
}

// 定义函数：矩阵相加相减
void MatSum(NewMat *tmp, NewMat *A, NewMat *B, int JudgeNum)
{
    if ((A->row != B->row) || (A->col != B->col))
    {
        // ! 这里需要一处报错，报错内容：(error:Not Match Matrix Dims)
        printf("error:矩阵维度错误！！！\n");
    }
    else
    {
        // IniNewMat(tmp, A->row, A->col);
        tmp->row = A->row;
        tmp->col = A->col;
        if (JudgeNum > 0)
        {
            for (int ind_row = 0; ind_row < A->row; ind_row++)
            {
                for (int ind_col = 0; ind_col < A->col; ind_col++)
                {
                    tmp->mat[ind_row][ind_col] = A->mat[ind_row][ind_col] + B->mat[ind_row][ind_col];
                }
            }
        }
        else if (JudgeNum < 0)
        {
            for (int ind_row = 0; ind_row < A->row; ind_row++)
            {
                for (int ind_col = 0; ind_col < A->col; ind_col++)
                {
                    tmp->mat[ind_row][ind_col] = A->mat[ind_row][ind_col] - B->mat[ind_row][ind_col];
                }
            }
        }
    }
}

// 定义函数：矩阵相乘
void MatMulti(NewMat *tmp, NewMat *A, NewMat *B)
{
    if (A->col != B->row)
    {
        // ! 这里需要一处报错，报错内容：(error:Not Match Matrix Dims)
        printf("矩阵相乘error:矩阵维度错误！！！\n");
    }
    else
    {
        // IniNewMat(tmp, A->row, B->col);
        // tmp->row = A->row;
        // tmp->col = B->col;
        for (int cal_row = 0; cal_row < A->row; cal_row++)
        {
            for (int cal_col = 0; cal_col < B->col; cal_col++)
            {
                double CalRes = 0.0;
                for (int cal_mid = 0; cal_mid < A->col; cal_mid++)
                {
                    CalRes += (A->mat[cal_row][cal_mid] * B->mat[cal_mid][cal_col]);
                }
                tmp->mat[cal_row][cal_col] = CalRes;
            }
        }
    }
}

// 定义函数：矩阵转置
void MatTranspose(NewMat *tmp, NewMat *A)
{
    // IniNewMat(tmp, A->col, A->row);
    // tmp->row = A->row;
    // tmp->col = A->col;
    for (int ind_i = 0; ind_i < A->row; ind_i++)
    {
        for (int ind_j = 0; ind_j < A->col; ind_j++)
        {
            tmp->mat[ind_j][ind_i] = A->mat[ind_i][ind_j];
        }
    }
}

// 矩阵求行列式
double getA(NewMat *Mat, int n)
{
    if (n == 1)
    {
        return Mat->mat[0][0];
    }
    double ans = 0;
    // double temp[Mat->row][Mat->col];
    NewMat temp;
    IniNewMat(&temp, Mat->row, Mat->col);
    // int j, j, k;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n - 1; j++)
        {
            for (int k = 0; k < n - 1; k++)
            {
                temp.mat[j][k] = Mat->mat[j + 1][(k >= i) ? k + 1 : k];
            }
        }
        double t = getA(&temp, n - 1);
        if (i % 2 == 0)
        {
            ans += Mat->mat[0][i] * t;
        }
        else
        {
            ans -= Mat->mat[0][i] * t;
        }
    }
    FreeMat(&temp);
    return ans;
}

// 矩阵求伴随矩阵
void getAStart(NewMat *Mat, int n, NewMat *AStar)
{
    if (n == 1)
    {
        Mat->mat[0][0] = 1;
        return;
    }
    // int i, j, k, t;
    NewMat temp;
    IniNewMat(&temp, Mat->row, Mat->col);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < n - 1; k++)
            {
                for (int t = 0; t < n - 1; t++)
                {
                    temp.mat[k][t] = Mat->mat[k >= i ? k + 1 : k][t >= j ? t + 1 : t];
                }
            }
            AStar->mat[j][i] = getA(&temp, n - 1);
            if ((i + j) % 2 == 1)
            {
                AStar->mat[j][i] = -AStar->mat[i][j];
            }
        }
    }
    FreeMat(&temp);
}

void getAInverse(NewMat *Res, NewMat *Mat)
{
    Res->row = Mat->row;
    Res->col = Mat->col;
    NewMat AStar;
    IniNewMat(&AStar, Mat->row, Mat->col);
    double a = getA(Mat, Mat->row);
    if (a == 0)
    {
        printf("矩阵求逆error:矩阵维度错误！！！\n");
        // printf("can not transform!\n");
    }
    else
    {
        getAStart(Mat, Mat->row, &AStar);
        for (int i = 0; i < Mat->row; i++)
        {
            for (int j = 0; j < Mat->row; j++)
            {
                Res->mat[i][j] = AStar.mat[i][j] / a;
            }
        }
    }
    FreeMat(&AStar);
}

void MatCopy(NewMat *CopyMat, NewMat *OriMat)
{
    // if ((CopyMat->row != OriMat->row) && (CopyMat->col != OriMat->col))
    // {
    //     printf("矩阵复制error:矩阵维度错误!!!\n");
    // }
    // else
    // {
    // IniNewMat(CopyMat, OriMat->row, OriMat->col);
    for (int i = 0; i < OriMat->row; i++)
    {
        for (int j = 0; j < OriMat->col; j++)
        {
            CopyMat->mat[i][j] = OriMat->mat[i][j];
        }
    }
    // }
}

double gaussianRand(double mean, double variance)
{
    static double U, V;
    static int phase = 0;
    double z;

    if (phase == 0)
    {
        U = rand() / (RAND_MAX + 1.0);
        V = rand() / (RAND_MAX + 1.0);
        z = sqrt(-1.0 * log(U)) * sin(2.0 * PI * V);
    }
    else
    {
        z = sqrt(-1.0 * log(U)) * cos(2.0 * PI * V);
    }
    phase = 1 - phase;
    z = mean + (z * pow(variance, 1 / 2));
    return z;
}

void IniNewMatRec(NewMatRec *IniMat, int rowNew, int colNew, int pageNew)
{
    IniMat->row = rowNew;
    IniMat->col = colNew;
    IniMat->page = pageNew;
    IniMat->mat = (double ***)malloc(rowNew * sizeof(double **));
    for (int ind_row = 0; ind_row < IniMat->row; ind_row++)
    {
        IniMat->mat[ind_row] = (double **)malloc(colNew * sizeof(double *));
        for (int ind_col = 0; ind_col < IniMat->col; ind_col++)
        {
            IniMat->mat[ind_row][ind_col] = (double *)malloc(pageNew * sizeof(double));
            for (int ind_page = 0; ind_page < IniMat->page; ind_page++)
            {
                IniMat->mat[ind_row][ind_col][ind_page] = 0.0;
            }
        }
    }
}

void FreeMatRec(NewMatRec *ResMat)
{
    for (int ind_row = 0; ind_row < ResMat->row; ind_row++)
    {
        for (int ind_col = 0; ind_col < ResMat->col; ind_col++)
        {
            free(ResMat->mat[ind_row][ind_col]);
        }
        free(ResMat->mat[ind_row]);
    }
    free(ResMat->mat);
}

// void FreeMat(NewMat *ResMat)
// {
//     // free(ResMat->row);
//     // free(ResMat->col);
//     for (int i = 0; i < ResMat->row; i++)
//         free(ResMat->mat[i]);
//     // puts("free Memory Done.");
//     free(ResMat->mat);
// }