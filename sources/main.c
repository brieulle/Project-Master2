#include"../entetes/lib.h"



int main(int argc, char** argv)
{

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
//
	int upperbound = 270000;	
	base primebase;
	preS_set F;
	mpz_t N, k, res;
	mpz_inits(N, k, res, NULL);
	mpz_set_ui(k,1);
	mpz_set_str(N, "18446744073709551617", 10);
	S_set S;

	cfrac frac;
	
	int i = factorisation(N, k, 270000, upperbound, primes);
	if(i == -1)
		printf("Try again :(\n");
	

	frac = expand(N, 1000, k);
//cfrac_display(frac);
	mpz_init_set_ui(k,1);
	primebase = gen_primeBase(upperbound, N, k, primes);
printf("La base est de taille %u\n",primebase.taille);
base_display(primebase);
	F = pairingSquare(frac.Q, frac.borne, primebase, N);
	S = S_set_affect(F, frac.A, N); 
	printf("Taille du S-set : %d\n", S.taille);
	isFactored(N, S, res);
	gmp_printf("Attention susPANSE : %Zd\n", res);
/////
printf("par ");
for(int i = 0; i < primebase.taille; i++)
	printf("%u ", primebase.tab[i]);
printf("\n");
for(int j = 0; j < F.nb_Q; j++)
{
for(int i = primebase.taille; i > -1; i--)
	gmp_printf("%u  ", mpz_tstbit(F.valQ[j], i));
printf("\n");
gmp_printf("Position du premier 1 : %u", mpz_scan1(F.valQ[j], 0));
printf("\n");
}

for(int i = 0; i < F.nb_Q; i++)
{
printf("h_%d  = ", i);
	for(int j = F.nb_Q; j > 0; j--)
		gmp_printf("%u ", mpz_tstbit(F.h[i], j));
printf("\n");
}
*/
	fclose(primes);
	return 0;
}
