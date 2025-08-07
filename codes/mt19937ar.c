/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/
#include <stdlib.h>
#include <math.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

typedef struct
{
	unsigned long *mt; /* the array for the state vector  */
	int mti;
	int snormal_state;
 	double snormal_next_value;
} MT_state;

MT_state *init_mt_state(void)
{
	MT_state *mt_state;
	mt_state = (MT_state*)malloc(sizeof(MT_state));
	if (mt_state == NULL) {
		printf("Memory error\n");
		exit(1);
	}
	mt_state->mt = (unsigned long*)malloc(N * sizeof(unsigned long));
	if (mt_state->mt == NULL) {
		printf("Memory error\n");
		exit(1);
	}
	mt_state->mti = N+1; /* mti==N+1 means mt[N] is not initialized */
	mt_state->snormal_state = 0;
	return mt_state;
}

void dealloc_mt_state(MT_state *mt_state)
{
	free(mt_state->mt);
	free(mt_state);
}

/* initializes mt[N] with a seed */
void init_genrand(MT_state *mt_state, unsigned long s)
{
	unsigned long *mt;
	int mti;
	
	mt = mt_state->mt;
	mti = mt_state->mti;
	
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
	
	mt_state->mti = mti;
	mt_state->snormal_state = 0;
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(MT_state *mt_state, unsigned long init_key[], int key_length)
{
    int i, j, k;
	unsigned long *mt;
	
	mt = mt_state->mt;
	
    init_genrand(mt_state, 19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(MT_state *mt_state)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
	unsigned long *mt;
	int mti;
	
	mt = mt_state->mt;
	mti = mt_state->mti;

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1) {   /* if init_genrand() has not been called, */
            init_genrand(mt_state, 5489UL); /* a default initial seed is used */
			mti = mt_state->mti;
		}

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
	
	mt_state->mti = mti;
	
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(MT_state *mt_state)
{
    return (long)(genrand_int32(mt_state)>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(MT_state *mt_state)
{
    return genrand_int32(mt_state)*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(MT_state *mt_state)
{
    return genrand_int32(mt_state)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(MT_state *mt_state)
{
    return (((double)genrand_int32(mt_state)) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(MT_state *mt_state) 
{ 
    unsigned long a=genrand_int32(mt_state)>>5, b=genrand_int32(mt_state)>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

double snormal_rand(MT_state *mt_state)
{
	double u1, u2;
	double z1, z2;
	double s;
	
	if(mt_state->snormal_state == 0){
		u1 = genrand_real3(mt_state);
		u2 = genrand_real3(mt_state);
		
		s = sqrt(-2.0*log(u1));
		
		z1 = s*cos(2.0*3.1415926535897932*u2);
		z2 = s*sin(2.0*3.1415926535897932*u2);
		
		mt_state->snormal_next_value = z2;
 		mt_state->snormal_state = 1;
 		return z1;
	}else{
		mt_state->snormal_state = 0;
		return mt_state->snormal_next_value;
	}
}