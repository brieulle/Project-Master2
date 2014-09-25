#include"../entetes/lib.h"

void cfrac_display(const cfrac frac)
{
	mpz_t tempz;
	mpz_init(tempz);

	printf("Valeurs de Q_n\n");
	for(int i = 0; i < frac.borne; i++)
	{
		printf("n = %d, Q_%d = ", i-1, i-1);
		gmp_printf("%Zd\n", frac.Q[i]);
	}

	printf("\nValeurs de g + P_n\n");
	for(int i = 0; i < frac.borne; i++)
	{
		printf("n = %d, g + P_%d = ", i, i);
		mpz_add(tempz, frac.g, frac.P[i]);
		gmp_printf("%Zd\n", tempz);
		mpz_set_ui(tempz, 0);
	}
		
	printf("\nValeurs de q_n\n");
	for(int i = 0; i < frac.borne; i++)
	{
		printf("n = %d, q_%d = ", i, i);
		gmp_printf("%Zd\n", frac.q[i]);
	}	


	printf("\nValeurs de r_n\n");
	for(int i = 0; i < frac.borne; i++)
	{
		printf("n = %d, r_%d = ", i-1, i-1);
		gmp_printf("%Zd\n", frac.r[i]);
	}	

	printf("\nValeurs de A_n\n");
	for(int i = 0; i < frac.borne; i++)
	{
		printf("n = %d, A_%d = ", i, i-1);
		gmp_printf("%Zd\n", frac.A[i]);
	}	
}


base gen_primeBase(const unsigned int UpperBound, const mpz_t N, const mpz_t k, FILE* primes)
{
	base primebase;
	int compteur, taille, temp_int;
	unsigned int* resultat = NULL;		//Contient les premiers qui formeront la base, critère : être un résidu quadratique modulo kN
	unsigned int* temp_res = NULL;			//Contient tous les nombres premiers strictement inférieur à UpperBound
	char* temp_char = NULL;			//Contient les char lors de la lecture du fichier
	mpz_t temp_mpz, temp_kN;

//	mpz_inits(temp_mpz, temp_kN, NULL);
	mpz_init(temp_mpz);
	mpz_init(temp_kN);
	taille = 0;	
	
	temp_char = (char*)malloc(10*sizeof(char));	//Le but étant de prendre des nombres premiers petits, il parait raisonnable de ne pas considérer des nombres premiers plus grand que 1000000, d'où le 6
	temp_res = (unsigned int*)malloc(UpperBound*sizeof(unsigned int));


	//On commence par récupérer les nombres premiers plus petit que UpperBound dans un fichier
	while(1)
	{
		fgets(temp_char, 10, primes);
		temp_res[taille] = atoi(temp_char);	
		taille++;

		if(temp_res[taille-1] >= UpperBound)	//Dès qu'on dépasse la borne, on s'arrête
			break;
	}

	temp_res = realloc(temp_res,(taille-1)*sizeof(unsigned int));

	resultat = (unsigned int*)malloc((taille-1)*sizeof(unsigned int));	//Au pire tous les premiers de tempr fonctionnent; à priori il n'y en aura que la moitié, environ
	
	mpz_set(temp_kN, k);
	mpz_mul(temp_kN, temp_kN, N);
	resultat[0] = 2;				//Nécessaire, la fonction mpz_legendre ne prend pas en compte le cas p = 2 et puis que tous les nombres impairs sont des carrés et que si les Q_n sont pairs ce serait bête de perdre leur division par deux et bien on l'ajoute de base !
	compteur = 1;

	for(int i = 1; i < taille; i++)	 
	{ 
		mpz_set_ui(temp_mpz, temp_res[i]);
		temp_int = mpz_legendre(temp_kN, temp_mpz);
		
		//Si le symbole de Legendre est nul, alors le nombre premier divise kN, on a alors trouvé un facteur
		if(temp_int == 0)
		{
			primebase.tab = (unsigned int*)malloc(sizeof(unsigned int));
			primebase.tab[0] = temp_res[i];
			primebase.taille = -1;

			return primebase;
		}
	
		//On ne retient que les premiers dont kN est un résidu quadratique en accord avec le théorème cité 
		if(temp_int == 1)	
		{
			resultat[compteur] = temp_res[i];
			compteur++;
		}
	}
	
	primebase.tab = (unsigned int*)realloc(resultat, compteur*sizeof(unsigned int));
	primebase.taille = compteur;

	return primebase;		//Penser à tester si la base est nulle ou non (au cas où le realloc foire)

}

void base_display(const base primebase)
{
	if(primebase.taille == -1)
		printf("félicitation, t'as trouvé un facteur ! : %u\n", primebase.tab[0]);
	else{
		for(int i = 0; i < primebase.taille; i++)
			printf("%u ", primebase.tab[i]);

		printf("\n taille = %u\n", primebase.taille);
	}
}






