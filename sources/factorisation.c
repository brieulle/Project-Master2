#include"../entetes/lib.h"


int factorisation(const mpz_t N, const mpz_t k, const long long unsigned int rang, const long long unsigned int UpperBound, mpz_t facteur)	//Fonction qui permet d'exhiber un facteur de N
{
	base primebase;		//La base de factorisation des Q_n	
	S_set S;			//La structure qui contiendra les paires A - Q
	preS_set preS;		//Étape intermédiaire
	cfrac frac;			//Le développement en fraction continue de N

/*	if(mpz_divisible_ui_p(N, 2)) //Inutile si on respect les données du problème
	{
		mpz_set_ui(facteur, 2);
		gmp_printf("Un facteur de %Zd est %d\n", N, 2);
		return 0;
	}*/

	primebase = gen_primeBase(UpperBound, N, k);
	if(primebase.taille == -1)			
	{
printf("m1\n");
		mpz_init_set(facteur, primebase.tab[0]);
		gmp_printf("Un facteur de %Zd est %Zd\n", N, facteur);
		return 0;
	}

	if(primebase.taille > rang)	//failsafe pour éviter les segfault si jamais le rang jusqu'où on développe la fraction continue est trop petit par rapport à la taille de la base
	{
		printf("Augmenter le rang de la fraction continue, au moins %llu\n", primebase.taille);
		return -3;
	}

	frac = expand(N, rang, k);
	if(frac.rang == 0)					//On vérifie qu'il n'y a pas eu d'erreur pendant le calcul des Q_n qui doivent toujours être < 2*sqrt(kN)
		return 1;
	preS = pairingSquare(frac, primebase, N);
	S = S_set_affect(preS, frac.A, N); 
	
	if(S.taille == 0)
		return -2;		

	isFactored(N, S, facteur);

	if(mpz_cmp_si(facteur, -1) == 0)
		return -1;
	else{
		gmp_printf("Un facteur de %Zd est %Zd\n", N, facteur);
		return 0;
	}
}


void boucleFactorisation(const mpz_t N, const mpz_t k, const long long unsigned int upperbound, const long long unsigned int rang, const long long unsigned int uplevel)
{
	int success = 0;
	long long unsigned int temp_rang, temp_upperbound;
	int factorisation_result = 0;
	mpz_t factor1, factor2, tempk;		

	mpz_inits(factor1, factor2, tempk, NULL);

	temp_rang = rang;
	temp_upperbound = upperbound;
	mpz_set(tempk, k);
//int temp = 0;
	while(success == 0 /*&& temp == 0*/)
	{
		factorisation_result = factorisation(N, tempk, temp_rang, temp_upperbound, factor1);

		switch(factorisation_result)
		{
			case -3:
				temp_rang = temp_upperbound;
				break;
			case -2:
				printf("Pas de S-set, changez k si vous gardez la même borne pour la base de factorisation ou vous pouvez augementer le rang de la fraction continue (ce qui est fait automatiquement).\n");
				temp_rang += uplevel;
				temp_upperbound += uplevel;
				break;
			case -1:
				printf("Pas de bon pgcd\n");
				temp_rang += uplevel;
				temp_upperbound += uplevel;
				break;
			case 0:
				success++;
				break;
			default:
				printf("Un problème s'est produit au niveau du calcul des Q_n.\n");
				return;
		}
//temp++;
	}

	if(success != 0)
	{
		mpz_divexact(factor2, N, factor1);
		gmp_printf("On a alors %Zd = %Zd.%Zd\n", N, factor1, factor2);
		mpz_clears(factor1, factor2, tempk, NULL);

		return;
	}else{
		mpz_clears(factor1, factor2, tempk, NULL);

		return;
	}
}

