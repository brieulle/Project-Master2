#include"../entetes/lib.h"

S_set S_set_affect(const preS_set F, const mpz_t* A, const unsigned int card_base)
{
	S_set S;
	unsigned int taille = 0;
	int* indice = (int*)malloc(F.nb_Q*sizeof(int)); //Tableau pour stocker l'indice des vecteurs qui fourniront un carré, il y en a AU MIEUX le nombre de Q_n qu'on a su factoriser; Essaie de voir si c'est vraiment utile, parce que d'ici ça ressemble à une idée abandonnée qui sert plus à grand chose...:w

	for(int i = 0; i < F.nb_Q; i++)	//On cherche les vecteurs mpz_t nuls, i.e. les produits des Q_n B-friable qui forment un carré
	{
		if(mpz_cmp_ui(F.valQ[i],0) == 0)
		{	
			indice[taille] = i;
			taille++;
		}
	}	
			
	if(taille == 0)
	{
		S.taille = 0;
		S.Qs = NULL;
		S.As = NULL;
	
		return S;
	}else{
		S.taille = taille;

		for(int i = 0; i < taille; i++)
		{
			mpz_inits(S.Qs[i], S.As[i],	NULL);

			for(int j = 0; j < card_base; j++)	//On multiplie les (-1)^jQ_j et les A_j selon la valeur du j-ième bit
			{
				switch(j)
				{
					case 0:
						if(mpz_tstbit(F.h[i],j) == 1)
						{
							//TODO : calculer les Q_i en utilisant l'algo exposé dans le papier
							mpz_set(S.As[i], A[F.n[j]]);	//Pense à revoir es indices, mais à priori c'est bon puisqu'il y a u déjà un décalage de -1 avec es A_n
						}else{
							mpz_set_ui(S.Qs[i],1);			//On initialise la valeur pour éviter le produit nul plus tard
							mpz_set_ui(S.As[i],1);
						}

						break;		
					default:
						if(mpz_tstbit(F.h[i],j) == 1)
						{
							mpz_mul(S.Qs[i], S.Qs[i], F.nQ[j]);
							mpz_mul(S.As[i], S.As[i], A[F.n[j]]);
						}
						
						break;
				}
				
			}

			if(mpz_tstbit(F.h[i],card_base) == 1)			//Si le bit de parité est à 1
				mpz_submul_ui(S.Qs[i] ,S.Qs[i],1);	//Alors on multiplie par -1
		}
	}

	return S;
}
