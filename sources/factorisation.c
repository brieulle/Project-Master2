#include"../entetes/lib.h"


int factorisation(const mpz_t N, const mpz_t k, unsigned int borne, unsigned int UpperBound, FILE* file)	//Fonction qui permet d'exhiber un facteur de N
{
	base primebase;		//La base de factorisation des Q_n	
	S_set S;			//La structure qui contiendra les paires A - Q
	preS_set preS;		//Étape intermédiaire
	cfrac frac;			//Le développement en fraction continue de N
	mpz_t facteur; 		//Le facteur de N
	
	mpz_init(facteur);

	if(mpz_divisible_ui_p(N, 2))
	{
		gmp_printf("Un facteur de %Zd est %d\n", N, 2);
		return 0;
	}

	primebase = gen_primeBase(UpperBound, N, k, file);
	if(primebase.taille == -1)			
	{
printf("m1\n");
		mpz_init_set_ui(facteur, primebase.tab[0]);
		gmp_printf("Un facteur de %Zd est %Zd\n", N, facteur);
		return 0;
	}

	frac = expand(N, borne, k);
	preS = pairingSquare(frac.Q, borne, primebase, N);
	S = S_set_affect(preS, frac.A, N); 
	
	if(S.taille == 0)
		return -1;		

	isFactored(N, S, facteur);

	if(mpz_cmp_si(facteur, -1) == 0)
		return -1;
	else{
		gmp_printf("Un facteur de %Zd est %Zd\n", N, facteur);
		return 0;
	}
}


void boucleFactorisation(const mpz_t N, const mpz_t k)
{
}	
