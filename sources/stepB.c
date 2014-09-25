#include"../entetes/lib.h"


preS_set pairingsquare(const mpz_t* A, const mpz_t* Q, const unsigned int borne, const base primebase, const mpz_t N)
{
	preS_set F;
	int j, k, pos;

	//On commence par factoriser les Q_n selon la base de factorisation, en éliminant les Q_n qui ne sont pas factorisable
	F = factQ(Q, borne, primebase, N);	//Retourne les Q_n B-friable

	//Début de la réduction de Gauss par la droite, le but est de trouver un produit tel que toutes les valuations soient paires

	j = primebase.taille+1;

	while(j >= 0)
	{
		k = 0;
		pos = 0;

		//Le but est de trouver le premier vecteur tel que la j ème composante (en partant de 0) soit égal à 1 <=> au pos-ième bit du mpz_t égal à 1 
		while((mpz_scan1(F.valQ[k],primebase.taille+1 - j) > pos) || k < F.nb_Q)
		{
			k++;
		}

		if(k >= F.nb_Q)	//Si aucun vecteur ne statisfait la condition on passe à la composante suivante et on change la position du bit
		{
			j--;
			pos++;
		}else{											//Sinon, on applique la réduction
			for(int i = k+1; i < F.nb_Q; i++)
			{
				if(mpz_tstbit(F.valQ[i], pos) == 1)		//Uniquement sur les vecteurs qui ont un 1 à la j-ième composante ou pos-ième position
				{
					mpz_xor(F.valQ[i], F.valQ[i], F.valQ[k]);
					mpz_xor(F.h[i], F.h[i], F.h[k]);
					j--;
					pos++;
				}else{
					j--;
					pos++;					
				}
			}	
		}		
	}	
	
	return F;	//On renvoie F dans le but de calculer le S-set dans l'étape C
}



preS_set factQ(const mpz_t* Q, const unsigned int borne, const base primebase, const mpz_t N)
{
	mpz_t* valQ = NULL;		//Tableau qui contiendra les valuations p-adique modulo 2 des Q_n factorisables
	mpz_t* h = NULL;		//Vecteur historique associé à chaque Q_n factorisables
	mpz_t* nQ = NULL;		//Les Q_n factorisables
	int* n = NULL;			//Tableau qui contient les Q_n factorisables
	int* truval = NULL;		//Tableau qui contiendra les vrai valuations des p (pour vérifier si Q_n est factorisé ou non)
	unsigned int nb_Q = 0;	//Le nombre de Q_n Factorisables
	int compteur, temp_int;
	mpz_t tempz;

	n = (int*)malloc(borne+2*sizeof(int));
	truval = (int*)malloc(primebase.taille*sizeof(int));
	valQ = (mpz_t*)malloc(borne+2*sizeof(mpz_t)); 
	h = (mpz_t*)malloc(borne+2*sizeof(mpz_t));
	nQ = (mpz_t*)malloc(borne+2*sizeof(mpz_t));
	compteur = 0;
	mpz_init(tempz);

	for(int i = 0; i < borne+2; i++)
	{
		mpz_init(valQ[compteur]);								//On initialise malgré tout, si jamais le Q_i n'est pas B-friable, on reprendra au même en droit

		for(int j = 0; j < primebase.taille; j++)					//Ce sera j+1 dans les vecteurs de F_2^r
		{
			truval[j] = 0;												//Valuation p-adique
			mpz_set(tempz, Q[i]);
		
			if(mpz_divisible_ui_p(tempz, primebase.tab[j]) != 0)	//Le but est de trouver (s'il y en a) les valuations p-adiques impaires
			{
				truval[j]++;
				mpz_divexact_ui(tempz, tempz, primebase.tab[j]);	
			
				while(mpz_divisible_ui_p(tempz, primebase.tab[j]) != 0)		//On compte le nombre de fois que Q_i est divisible par p
				{
					truval[j]++;
					mpz_divexact_ui(tempz, tempz, primebase.tab[j]);	
				}
			}

			if(truval[j]&1)
				mpz_setbit(valQ[compteur], primebase.taille - (j+1));							//Si la valuation est impaire, on met le r-1-ème bit à 1 <=> la j+1ème composante à 1	

		}
		
				




		}
	}	
			

}





 
