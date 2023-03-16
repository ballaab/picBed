#include "/usr/local/include/pbc/pbc.h"
#include "/usr/local/include/pbc/pbc_test.h"
#include <time.h>
#include <math.h>

#define M_len 100
#define C_len 100
#define PK_len 100
clock_t start, stop, start1, stop1, start2, stop2;
double duration, duration1, duration2;
int main(int argc, char **argv)
{

    pairing_t pairing;
    char param[2048];
    size_t count = fread(param, 1, 2048, stdin);
    if (!count)
        pbc_die("input error");
    pairing_init_set_buf(pairing, param, count);

    // Initilize
    element_t g, h;
    element_t msk, mpk;
    element_init_G1(g, pairing);
    element_init_G2(h, pairing);
    element_init_Zr(msk, pairing);
    element_init_G1(mpk, pairing);
    element_random(g);
    element_random(h);
    element_random(msk);
    element_pow_zn(mpk, g, msk);


    // KeyGen
    element_t sk, pk, *PK, temp1;
    PK = (element_t *)pbc_malloc(PK_len * sizeof(element_t));

    element_init_Zr(temp1, pairing);
    element_init_Zr(sk, pairing);
    element_init_G2(pk, pairing);
    element_random(sk);
    //temp1 = sk;
    element_set(sk, temp1);
    element_init_G2(PK[0], pairing);
    element_pow_zn(PK[0], h, sk);     //PK[0]=h^x
    for(int i = 1; i < 2 * PK_len; i++){
        element_init_G2(PK[i], pairing);
        element_mul(temp1, temp1, sk);
        element_pow_zn(PK[i], h, temp1);
    }

    element_t token1, temp2;
    element_init_G2(token1, pairing);
    element_init_G1(temp2, pairing);  
    element_from_hash(temp2, pk, 10);
    element_pow_zn(token1, temp2, msk);  //token1=h1(pk)^a  h1-> G2

    element_clear(temp1);
    element_clear(temp2);

    //TKeyGen
    element_t beta;
    element_init_Zr(beta, pairing);
    element_random(beta);

    // TagGen
    // 测试用：生成Zr M[100]，进行审计
    element_t *M, token2, *T, bunch_tag;
    element_t temp3, temp4, temp5, xk;

    element_init_Zr(temp2, pairing);
    element_init_G1(temp3, pairing);
    element_init_Zr(temp4, pairing);
    element_init_G1(temp5, pairing);
    element_init_G1(bunch_tag, pairing);
    element_init_Zr(xk, pairing);
    element_random(temp4);
    element_random(bunch_tag);

    M = (element_t *)pbc_malloc(M_len * sizeof(element_t));
    T = (element_t *)pbc_malloc(M_len * sizeof(element_t));

    //xk = sk;
    element_set(sk, xk);
    start = clock();
    for (int i = 0; i < M_len; i++)
    {
        element_init_Zr(M[i], pairing);
        element_init_G1(T[i], pairing);

        element_random(M[i]);
        element_pow_zn(temp3, g, M[i]); // temp3=g^mi
        element_from_hash(temp5, temp4, 1);  //temp5=H2(i)
        element_mul(temp3, temp3, temp5);   //temp3=(H2(i)*g^mi)
        element_pow_zn(T[i], temp3, xk);   //T[i]= temp3^xk
        element_mul(xk, xk, sk);
        element_mul(bunch_tag, bunch_tag, T[i]);
    }

    stop = clock();
    duration = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("%f seconds\n", duration / 10);

    element_clear(temp2);
    element_clear(temp3);
    element_clear(temp4);
    element_clear(temp5);

    //Challenge
    element_t *C;
    C = (element_t *)pbc_malloc(sizeof(element_t) * C_len);
    // b = pbc_malloc(sizeof(element_t) * chal_len);
    for (int i = 0; i < C_len; i++)
    {
        element_init_Zr(C[i], pairing);
        element_random(C[i]);
        // element_printf("a[i]=%B\n", a[i]);
    }

    printf("chal over\n");

    
    //Respond
    element_t temp6, temp7, temp8, temp9,temp10, temp11, temp12;
    element_t *EI, mu, sig;
    EI = (element_t *)pbc_malloc(M_len * sizeof(element_t));

    element_init_G1(temp6, pairing);
    element_init_Zr(temp7, pairing);
    element_init_G1(temp8, pairing);
    element_init_Zr(temp9, pairing);
    element_init_G2(temp10, pairing);
    element_init_Zr(temp11, pairing);
    element_init_GT(temp12, pairing);
    element_init_Zr(mu, pairing);
    element_random(mu);
    element_random(temp7);

    element_init_GT(sig, pairing);
    element_random(sig);

    stop1 = clock();
    //temp9 = sk;
    element_set(sk, temp9);
    for(int i = 0; i < M_len; i++){
        element_mul(temp9, temp9, sk);
    }
    element_mul(temp9, temp9, beta);  //temp9 = beta x^{n+1}
    element_pow_zn(temp10, h, temp9);  //temp10 = h^temp9

    for(int i = 0; i < M_len; i++){
        element_pow_zn(temp6, g, M[i]); // temp6=g^mi
        element_from_hash(temp8, temp7, 1);  //temp8=H2(i)
        element_mul(temp6, temp6, temp8);   //temp6=(H2(i)*g^mi)
        element_init_GT(EI[i], pairing);
        pairing_apply(EI[i], temp6, temp10, pairing);
    }

    for(int i = 0; i < C_len; i++){
        element_mul(temp11, C[i], M[i]);
        element_add(mu, mu, temp11);
        element_pow_zn(temp12, EI[i], C[i]);  
        element_mul(sig, temp12, sig);
    }
    duration1 = (double)(stop1 - start1) / CLOCKS_PER_SEC;
    printf("%f seconds\n", duration1 / 10);
    element_clear(temp6);
    element_clear(temp7);
    element_clear(temp8);
    element_clear(temp9);
    element_clear(temp10);
    element_clear(temp11);
    element_clear(temp12);

    //Verify
    element_t temp13, temp14, temp15, temp16,temp17, temp18, temp19;
    element_t validres, E;

    element_init_Zr(temp13, pairing);
    element_init_G2(temp14, pairing);
    element_init_Zr(temp15, pairing);
    element_init_G1(temp16, pairing);
    element_init_G1(temp17, pairing);
    element_init_G1(temp18, pairing);
    element_init_G1(temp19, pairing);
    element_init_GT(E, pairing);
    element_random(temp15);
    element_random(temp19);

    start2 = clock();
    element_from_hash(temp16, temp15, 1);  //temp16=H2(i)
    element_pow_zn(temp17, g, mu); //temp17=g^mu
    //temp13 = sk;
    element_set(sk, temp13);
    for(int i = 0; i < M_len; i++){
        element_mul(temp13, temp13, sk);
    }
    element_mul(temp13, temp13, beta);  //temp13 = beta x^{n+1}
    element_pow_zn(temp14, h, temp13);  //temp14 = h^temp9

    for(int i = 0; i < C_len; i++){
        element_mul(temp18, temp16, C[i]);
        element_mul(temp19, temp19, temp18);
    }
    element_mul(temp19, temp19, temp17);
    pairing_apply(E, temp19, temp14, pairing);

    //哈希验证

   
    stop2 = clock();
    duration2 = (double)(stop2 - start2) / CLOCKS_PER_SEC;
    printf("duration2: %f seconds\n", duration2 / 10);

    element_clear(temp13);
    element_clear(temp14);
    element_clear(temp15);
    element_clear(temp16);
    element_clear(temp17);
    element_clear(temp18);
    element_clear(temp19);
    element_clear(validres);

    element_clear(g);
    element_clear(msk);
    element_clear(mpk);
    element_clear(token1);
    element_clear(token2);
    element_clear(sk);
    element_clear(pk);
    element_clear(h);
    element_clear(mu);
    element_clear(sig);
    element_clear(beta);
    element_clear(bunch_tag);
    element_clear(xk);
    free(M);
    free(C);
    free(PK);
    free(T);
    free(EI);

}