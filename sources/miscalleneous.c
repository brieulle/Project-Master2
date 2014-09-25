#include"../entetes/lib.h"

void cfrac_display(const cfrac frac)
{
	mpz_t tempz;
	mpz_init(tempz);

	printf("Valeurs de Q_n\n");
	for(long long int i = 0; i < frac.rang; i++)
	{
		printf("n = %lld, Q_%lld = ", i-1, i-1);
		gmp_printf("%Zd\n", frac.Q[i]);
	}

	printf("\nValeurs de g + P_n\n");
	for(long long int i = 0; i < frac.rang; i++)
	{
		printf("n = %lld, g + P_%lld = ", i, i);
		mpz_add(tempz, frac.g, frac.P[i]);
		gmp_printf("%Zd\n", tempz);
		mpz_set_ui(tempz, 0);
	}
		
	printf("\nValeurs de q_n\n");
	for(long long int i = 0; i < frac.rang; i++)
	{
		printf("n = %lld, q_%lld = ", i, i);
		gmp_printf("%Zd\n", frac.q[i]);
	}	


	printf("\nValeurs de r_n\n");
	for(long long int i = 0; i < frac.rang; i++)
	{
		printf("n = %lld, r_%lld = ", i-1, i-1);
		gmp_printf("%Zd\n", frac.r[i]);
	}	

	printf("\nValeurs de A_n\n");
	for(long long int i = 0; i < frac.rang; i++)
	{
		printf("n = %lld, A_%lld = ", i, i-1);
		gmp_printf("%zd\n", frac.A[i]);
	}	
	
	mpz_clear(tempz);
}


base gen_primeBase(const long long unsigned int upperbound, const mpz_t N, const mpz_t k)
{
	FILE* primes = fopen("primelist.txt", "r");		//Fichier contenant les nombres premiers jusqu'à dix millions
	base primebase;
	long long int compteur, taille, temp_int;
	mpz_t* resultat = NULL;		//Contient les premiers qui formeront la base, critère : être un résidu quadratique modulo kn
	mpz_t* temp_res = NULL;			//Contient tous les nombres premiers strictement inférieur à upperbound
	char* temp_char = NULL;			//Contient les char lors de la lecture du fichier
	mpz_t temp_mpz, temp_kN;

	mpz_inits(temp_mpz, temp_kN, NULL);
	taille = 0;	
	
	temp_char = (char*)malloc(10*sizeof(char));	//Le but étant de prendre des nombres premiers petits, il parait raisonnable de ne pas considérer des nombres premiers plus grand que 1000000, d'où le 6
	temp_res = (mpz_t*)malloc(upperbound*sizeof(mpz_t));

	//On commence par récupérer les nombres premiers plus petit que upperbound dans un fichier
	while(1)
	{
		fgets(temp_char, 10, primes);
		mpz_init_set_si(temp_res[taille], atoll(temp_char));	
		taille++;

		if(mpz_cmp_ui(temp_res[taille-1], upperbound-1)>=0)	//Dès qu'on dépasse la borne, on s'arrête; upperbound-1 parce que je ne veux pas récupérer le premier égal à  upperbound s'il existe
			break;
	}

	fclose(primes);
	resultat = (mpz_t*)malloc((taille-1)*sizeof(mpz_t));	//Au pire tous les premiers de tempr fonctionnent; à priori il n'y en aura que la moitié, environ
	mpz_set(temp_kN, k);
	mpz_mul(temp_kN, temp_kN, N);
	mpz_init_set_ui(resultat[0], 2);				//Nécessaire, la fonction mpz_legendre ne prend pas en compte le cas p = 2 et puis que tous les nombres impairs sont des carrés et que si les q_n sont pairs ce serait bête de perdre leur division par deux et bien on l'ajoute de base !
	compteur = 1;

	for(long long int i = 1; i < taille; i++)	 
	{ 
		temp_int = mpz_legendre(temp_kN, temp_res[i]);
		
		//Si le symbole de Legendre est nul, alors le nombre premier divise kN, on a trouvé un facteur; il faut aussi que p ne divise pas k, sinon cela peut générer des erreurs (on trouverait alors un facteur de k)
		if(temp_int == 0 && mpz_divisible_p(k, temp_res[i]) == 0)
		{
			primebase.tab = (mpz_t*)malloc(sizeof(mpz_t));
			mpz_init_set(primebase.tab[0], temp_res[i]);
			primebase.taille = -1;

			return primebase;
		}
	
		//On ne retient que les premiers dont kN est un résidu quadratique en accord avec le théorème cité 
		if(temp_int == 1)	
		{
			mpz_init_set(resultat[compteur], temp_res[i]);
			compteur++;
		}
	}
	
	primebase.tab = resultat;
	primebase.taille = compteur;
	mpz_clears(temp_mpz, temp_kN, NULL);

	return primebase;		//Penser à tester si la base est nulle ou non (au cas où le realloc foire)

}

void base_display(const base primebase)
{
	if(primebase.taille == -1)
		gmp_printf("Félicitation, t'as trouvé un facteur ! : %Zd\n", primebase.tab[0]);
	else{
		for(long long int i = 0; i < primebase.taille; i++)
			gmp_printf("%Zd ", primebase.tab[i]);

		printf("\n taille = %llu\n", primebase.taille);
	}
}

