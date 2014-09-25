#include"../entetes/lib.h"


preS_set pairingSquare(cfrac frac, const base primebase, const mpz_t N)
{
	preS_set F;
	long long unsigned int k, pos;

	//On commence par factoriser les Q_n selon la base de factorisation, en éliminant les Q_n qui ne sont pas factorisable
	F = factQ(frac.Q, frac.rang, primebase, N);	//Retourne les Q_n B-friable

	//Début de la réduction de Gauss par la droite, le but est de trouver un produit tel que toutes les valuations soient paires

	pos = 0;

	while(pos < primebase.taille)
	{
		k = 0;

		//Le but est de trouver le premier vecteur tel que la j ème composante (en partant de 0) soit égal à 1 <=> au pos-ième bit du mpz_t égal à 1 
		for(long long int i = 0; i < F.nb_Q; i++)
		{
			if(mpz_scan1(F.valQ[i], pos) == pos)
			{
				k=i;
				break;
			}
			k++;
		}	

		if(k >= F.nb_Q)	//Si aucun vecteur ne statisfait la condition on passe à la composante suivante et on change la position du bit
			pos++;
		else{											//Sinon, on applique la réduction
			for(long long int i = k+1; i < F.nb_Q; i++)
			{
				if(mpz_tstbit(F.valQ[i], pos) == 1)		//Uniquement sur les vecteurs qui ont un 1 à la j-ième composante ou pos-ième position
				{
					mpz_xor(F.valQ[i], F.valQ[i], F.valQ[k]);
					mpz_xor(F.h[i], F.h[i], F.h[k]);
				}
			}	
			pos++;	//on passe à la position suivante
		}		
	}

	return F;	//On renvoie F dans le but de calculer le S-set dans l'étape C
}

//Fonction dont le but est de factoriser les Q_n sur la base de factorisation générée précédemment
preS_set factQ(mpz_t* Q, const long long unsigned int rang, const base primebase, const mpz_t N)
{
	preS_set res;
	mpz_t* valQ = NULL;		//Tableau qui contiendra les valuations p-adique modulo 2 des Q_n factorisables
	mpz_t* h = NULL;		//Vecteur historique associé à chaque Q_n factorisables
	mpz_t* nQ = NULL;		//Les Q_n factorisables
	long long int* n = NULL;			//Tableau qui contient les indices Q_n factorisables parmis les Q_n de base
	long long int* truval = NULL;		//Tableau qui contiendra les vrai valuations p-adiques (pour vérifier si Q_n est factorisé ou non)
	long long int nb_Q, temp_int, compteur;
	mpz_t tempz, tempz_2;

	n = (long long int*)malloc((rang+2)*sizeof(long long int));
	truval = (long long int*)malloc(primebase.taille*sizeof(long long int));
	valQ = (mpz_t*)malloc((rang+2)*sizeof(mpz_t)); 
	h = (mpz_t*)malloc((rang+2)*sizeof(mpz_t));
	nQ = (mpz_t*)malloc((rang+2)*sizeof(mpz_t));
	nb_Q = 0;
	compteur = 0;
	mpz_inits(tempz, tempz_2, NULL);

	for(long long int i = 0; i < rang+2; i++)
	{
		mpz_init(valQ[compteur]);	//On initialise malgré tout, si jamais le Q_i n'est pas B-friable, on reprendra au même en droit

		for(long long int j = 0; j < primebase.taille; j++)					//Ce sera j+1 dans les vecteurs de F_2^r
		{
			truval[j] = 0;												//Valuation p-adique
			mpz_set(tempz, Q[i]);

			if(mpz_divisible_p(tempz, primebase.tab[j]) != 0)	//Le but est de trouver (s'il y en a) les valuations p-adiques impaires
			{
				truval[j]++;
				mpz_divexact(tempz, tempz, primebase.tab[j]);	
			
				while(mpz_divisible_p(tempz, primebase.tab[j]) != 0)		//On compte le nombre de fois que Q_i est divisible par p
				{
					truval[j]++;
					mpz_divexact(tempz, tempz, primebase.tab[j]);	
				}
			}
		}
		


		//On remultiplie les facteurs trouvés dans le but de voir si Q_i est complètement factorisable sur la base
		mpz_set_ui(tempz, 1);
		for(long long int j = 0; j < primebase.taille; j++)
		{			
			if(truval[j] > 0)					//Si la valuation p-adique est plus grande que 1
			{
				mpz_pow_ui(tempz_2, primebase.tab[j], truval[j]);
				mpz_mul(tempz, tempz, tempz_2);
			}
		}

		//Si Q_i est complètement factorisable sur la base alors on l'ajoute au nouveau tableau
		if(mpz_cmp(tempz, Q[i]) == 0 && mpz_cmp_ui(tempz, 1) != 0)
		{
			if((i-1)&1 && i!=0)								//Se rappeler que les indices sont décalés depuis l'étape A
				mpz_setbit(valQ[compteur],primebase.taille);
			for(long long int j = 0; j < primebase.taille; j++)
			{
				if(truval[j]&1)
					mpz_setbit(valQ[compteur], primebase.taille - j - 1);
			}	
			mpz_set(nQ[compteur], Q[i]);
			n[compteur] = i;
			compteur++;

		}
	}	
	
	if(compteur > 0)
	{
		for(long long int i = 0; i < compteur; i++)
		{
			mpz_init(h[i]);
			mpz_setbit(h[i], compteur - i);
		}
		
		res.valQ = realloc(valQ, compteur*sizeof(mpz_t));
		res.n = realloc(n, compteur*sizeof(long long int));
		res.h = realloc(h, compteur*sizeof(mpz_t));
		res.nQ = realloc(nQ, compteur*sizeof(mpz_t));
		res.nb_Q = compteur;

		mpz_clears(tempz, tempz_2, NULL);
		
		return res;
	}else{							//Si on ne trouve aucun Q_n factorisable
		res.nb_Q = 0;

		mpz_clears(tempz, tempz_2, NULL);
		
		return res;
	}
}
 
