#include"../entetes/lib.h"



int main(int argc, char** argv)
{
/*
	mpz_t N, k;
	mpz_inits(N, k, NULL);
	mpz_set_ui(k,1);
	mpz_set_str(N, "13290059", 10);

	cfrac frac;
	

	frac = expand(N, 4, k);
	cfrac_display(frac);

*/
	FILE* primes = fopen("Nombres premiers.txt", "r");
/*
	mpz_t z, y;
	mpz_init_set_ui(z, 7);
	mpz_init(y);
	int a,b,c,d;
	a = mpz_tstbit(z,0);
	b = mpz_tstbit(z,1);
	c = mpz_tstbit(z,2);
	d = mpz_tstbit(z,3);
	
	printf("%d%d%d%d\n",d,c,b,a);
	mp_bitcnt_t i = mpz_scan0(z, 1);
	gmp_printf("%u\n",i);
	gmp_printf("%Zd\n",y);
*/
	mpz_t N,k,a;
	int upperbound = 50;	
	base primebase;

	mpz_init_set_ui(k,1);
	mpz_init(N);
	mpz_set_str(N, "132901", 10);
	primebase = gen_primeBase(upperbound, N, k, primes);
	base_display(primebase);

	fclose(primes);
	return 0;
}
