#include"../entetes/lib.h"

S_set S_set_affect(const preS_set F, mpz_t* A, const mpz_t N)
{
	S_set S;
	unsigned int taille_0, taille_Q;				//taille_0 = le nombre de h_i = 0, taille_Q = le nombre de Q_j "dans" les h_j
	int* indice = (int*)malloc(F.nb_Q*sizeof(int)); //Tableau pour stocker l'indice des vecteurs qui fourniront un carré, il y en a AU MIEUX le nombre de Q_n qu'on a su factoriser; Essaie de voir si c'est vraiment utile, parce que d'ici ça ressemble à une idée abandonnée qui sert plus à grand chose...:w
	mpz_t R, X, Q;	//Les R, X et Q tel que décrit dans la méthode de racine carré mod N
	mpz_t* tempQ = (mpz_t*)malloc(F.nb_Q*sizeof(mpz_t));	//Contiendra les Q_j qui forment les carrés des paires A - Q; ils seront déterminés selon si leur indice apparait dans les h_i = 0
	mpz_t* Qs = NULL;
	mpz_t* As = NULL;
	mpz_inits(R, X, Q, NULL);
	taille_0 = 0;

	for(int i = 0; i < F.nb_Q; i++)	//On cherche les vecteurs mpz_t nuls, i.e. les produits des Q_n B-friable qui forment un carré
	{
		if(mpz_cmp_ui(F.valQ[i],0) == 0)
		{	
			indice[taille_0] = i;
			taille_0++;
		}
	}	
			
	if(taille_0 == 0)
	{
		S.taille = 0;
		Qs = NULL;
		S.As = NULL;
	
		return S;
	}else{
		Qs = (mpz_t*)malloc(taille_0*sizeof(mpz_t));
		As = (mpz_t*)malloc(taille_0*sizeof(mpz_t));

		for(int i = 0; i < taille_0; i++)
		{
			mpz_inits(Qs[i], As[i],	NULL);

			//On commence par le calcul de As qui est assez direct 
			for(int j = 0; j < F.nb_Q; j++)	
			{
				switch(j)
				{
					case 0:
						if(mpz_tstbit(F.h[indice[i]],j) == 1)		//Le but est de rajouter les A_i-1 correspondant au Q_i selon les valeurs des bit de h[i] associé aux valQ = 0
						{
							mpz_init_set(As[i], A[F.n[j]]);
						}else{
							mpz_init_set_ui(As[i],1);
						}

						break;		
					default:
						if(mpz_tstbit(F.h[indice[i]],j) == 1)
						{
							mpz_mul(As[i], As[i], A[F.n[j]]);
							mpz_mod(As[i], As[i], N);
						}
						
						break;
				}
				
			}

			//Calcul de Qs associé
			taille_Q = 0;

			for(int j = F.nb_Q; j > -1; j--)	//On commence à remplir le tableau des Q_j associés aux h_i
			{
				if(mpz_tstbit(F.h[indice[i]],j) == 1)
				{
					mpz_init_set(tempQ[taille_Q], F.nQ[j]);
					taille_Q++;
				}
			}
			
			//Début du Square Root Procedure
			//Étape 1
			mpz_set(R, tempQ[0]); 
			mpz_set_ui(Qs[i], 1);

			for(int j = 1; j < taille_Q; j++)
			{
				//Étape 2
				mpz_gcd(X, R, tempQ[j]);

				//Étape 3
				mpz_mul(Qs[i], Qs[i], X);
				mpz_mod(Qs[i], Qs[i], N);

				//Étape 4
				mpz_divexact(R, R, X);
				mpz_mul(R, R, tempQ[j]);
				mpz_divexact(R, R, X);
			}
	
			//Étape 7
			mpz_sqrt(X, R);
			
			//Étape 8
			mpz_mul(Qs[i], Qs[i], X);
			mpz_mod(Qs[i], Qs[i], N);
		}
		
		S.Qs = Qs;
		S.As = As;
		S.taille = taille_0;

		return S;
	}

}

void isFactored(const mpz_t N, S_set S, mpz_t res)	//Fonction qui renvoie le facteur si on en trouve un ou -1 si on en a pas trouvé
{
	mpz_t sqrA, sqrQ, deux, PGCD, temp;
	mpz_inits(sqrA, sqrQ, PGCD, temp, NULL);
	mpz_init_set_ui(deux, 2);
	int success = 0;			//Savoir si la factorisation a fonctionné

	for(int i = 0; i < S.taille; i++)
	{
		mpz_powm(sqrA, S.As[i], deux, N);
		mpz_powm(sqrQ, S.Qs[i], deux, N);

		if(mpz_cmp(sqrA, sqrQ) == 0)
			continue;					//Si les deux termes sont égaux, la méthode ne fonctionne pas, alors on passe à la prochaine paire
		else{
			mpz_sub(temp, sqrA, sqrQ);	
			mpz_mod(temp, temp, N);
			mpz_gcd(PGCD, N, temp);		//Si jamais la paire est bonne, le pgcd sera alors le facteur cherché
			
			if(mpz_cmp_ui(PGCD, 1)==0)	//Si le candidat est premier avec N, alors on a pas de facteur et on continue
				continue;
			else{				//On a trouve une paire qui marche alors on sort de la boucle
				success++;	
				break;
			}
		}
	}
		
	if(success == 0)
		mpz_set_si(res, -1);
	else
		mpz_set(res, PGCD);

}







