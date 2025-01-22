/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
  int i, j, i1, j1, blk = 8;
  int col,row;

  switch(M){
    case 32:
      //diagonal blocks, first move then to the block next to it and then copy to the correct block.
      for (i = 0; i < M; i += blk)
        if(i!=(M-blk)){
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1+blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1+blk]; 
        }else{
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1-blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1-blk]; 
        } 
      for (i = 0; i < M; i += blk)
        for (j = 0; j < M; j += blk){
          if(i!=j){
            for (i1 = i; i1 < i + blk; i1++)
              for (j1 = j; j1 < j + blk;  j1++)
                B[j1][i1] =  A[i1][j1]; 
          }    
        }
      break;
    case 64:
      blk = 8;
      //help from : "https://github.com/simnalamburt/snucse/blob/main/System%20Programming/cachelab/part_b/trans.c"
      //modified to reach hitrate 1292 (< 1300)
      for (i = 0; i < M; i += blk)
        if(i!=(M-blk)){
          // B의 윗쪽 반
          for (int k = 0; k < 8; ++k) {
            int *t = &A[i + k][i];
            int a = t[0], b = t[1], c = t[2], d = t[3];
            t = &B[i][i + k + blk];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 아랫쪽 반
          for (int k = 7; k >= 0; --k) {
            int *t = &A[i + k][i + 4];
            int a = t[0], b = t[1], c = t[2], d = t[3];
            t = &B[i][i + k + blk*2];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 윗쪽 반
          for (int k = 0; k < 8; ++k) {
            int *t =&B[i][i + k + blk];
            int a = t[0], b = t[64], c = t[128], d = t[192];
            t = &B[i][i + k];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 아랫쪽 반
          for (int k = 7; k >= 0; --k) {
            int *t = &B[i][i + k + blk*2];
            int a = t[0], b = t[64], c = t[128], d = t[192];
            t = &B[i + 4][i + k];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }     
        }else{
          // B의 윗쪽 반
          for (int k = 0; k < 8; ++k) {
            int *t = &A[i + k][i];
            int a = t[0], b = t[1], c = t[2], d = t[3];
            t = &B[i][i + k - blk];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 아랫쪽 반
          for (int k = 7; k >= 0; --k) {
            int *t = &A[i + k][i + 4];
            int a = t[0], b = t[1], c = t[2], d = t[3];
            t = &B[i][i + k - blk*2];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 윗쪽 반
          for (int k = 0; k < 8; ++k) {
            int *t =&B[i][i + k - blk];
            int a = t[0], b = t[64], c = t[128], d = t[192];
            t = &B[i][i + k];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

          // B의 아랫쪽 반
          for (int k = 7; k >= 0; --k) {
            int *t = &B[i][i + k - blk*2];
            int a = t[0], b = t[64], c = t[128], d = t[192];
            t = &B[i + 4][i + k];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }   
        }
      for(col = 0; col < N; col += blk)
      {
        for(row = 0; row < M; row += blk){
          if(row!=col){
            int *t = &A[col][row + 4];
            int a = t[0], b = t[1], c = t[2], d = t[3];
            // B의 윗쪽 반
            for (int k = 0; k < 8; ++k) {
              int *t = &A[col + k][row];
              int a = t[0], b = t[1], c = t[2], d = t[3];
              t = &B[row][col + k];
              t[0] = a; t[64] = b; t[128] = c; t[192] = d;
            }

            // B의 아랫쪽 반
            for (int k = 7; k > 0; --k) {
              int *t = &A[col + k][row + 4];
              int a = t[0], b = t[1], c = t[2], d = t[3];
              t = &B[row + 4][col + k];
              t[0] = a; t[64] = b; t[128] = c; t[192] = d;
            }
            t = &B[row + 4][col];
            t[0] = a; t[64] = b; t[128] = c; t[192] = d;
          }

        }
      }

      break;

    default:
      //tag
      blk = 17;  // Block size
      for (i = 0; i < N; i += blk) {
        for (j = 0; j < M; j += blk) {
          // Transpose blk x blk sub-matrix
          for (i1 = i; i1 < (i + blk) && i1 < N; i1++) {
            for (j1 = j; j1 < (j + blk) && j1 < M; j1++) {
              B[j1][i1] = A[i1][j1];  // Transpose element

            }
          }
        }
      }
      break;
  }

}

char transpose_submit2_desc[] = "Transpose submission2";
void transpose_submit2(int M, int N, int A[N][M], int B[M][N])
{
  int i, j, i1, j1, blk;
  switch(M){
    case 32:
      blk = 8;
      //diagonal blocks, first move then to the block next to it and then copy to the correct block.
      for (i = 0; i < M; i += blk)
        if(i!=(M-blk)){
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1+blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1+blk]; 
        }else{
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1-blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1-blk]; 
        } 
      //all the blocks left
      for (i = 0; i < M; i += blk)
        for (j = 0; j < M; j += blk){
          if(i!=j){
            for (i1 = i; i1 < i + blk; i1++)
              for (j1 = j; j1 < j + blk;  j1++)
                B[j1][i1] =  A[i1][j1]; 
          }    
        }
      break;
    case 64:
      blk = 4;
      //diagonal blocks, first move then to the block next to it and then copy to the correct block.
      for (i = 0; i < M; i += blk)
        if(i!=(M-blk)){
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1+blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1+blk]; 
        }else{
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1-blk] =  A[i1][j1]; 
          for (i1 = i; i1 < i + blk; i1++)
            for (j1 = i; j1 < i + blk;  j1++)
              B[j1][i1] = B[j1][i1-blk]; 
        } 
      for (i = 0; i < M; i += blk)
        for (j = 0; j < M; j += blk){
          if(i!=j){
            for (i1 = i; i1 < i + blk; i1++)
              for (j1 = j; j1 < j + blk;  j1++)
                B[j1][i1] =  A[i1][j1]; 
          }    
        }
      break;

    default:
      //tag
      blk = 17;  // Block size
      for (i = 0; i < N; i += blk) {
        for (j = 0; j < M; j += blk) {
          // Transpose blk x blk sub-matrix
          for (i1 = i; i1 < (i + blk) && i1 < N; i1++) {
            for (j1 = j; j1 < (j + blk) && j1 < M; j1++) {
              B[j1][i1] = A[i1][j1];  // Transpose element

            }
          }
        }
      }
      break;
  }

}


/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
  int i, j, tmp;

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; j++) {
      tmp = A[i][j];
      B[j][i] = tmp;
    }
  }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
  /* Register your solution function */
  registerTransFunction(transpose_submit, transpose_submit_desc); 
  registerTransFunction(transpose_submit2, transpose_submit2_desc); 

  /* Register any additional transpose functions */
  registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < M; ++j) {
      if (A[i][j] != B[j][i]) {
        return 0;
      }
    }
  }
  return 1;
}

