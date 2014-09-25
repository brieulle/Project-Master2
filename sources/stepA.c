#include"../entetes/lib.h"

//Le but de la fonction est de calculer le développement en fractions continue jusqu'à un certain rang de racine carrée de kN et de stocker les couples (A_n-1,Q_n) comme décrit dans la section (à venir)
cfrac expand(const mpz_t N, const long long unsigned int rang, const mpz_t k)
{
	cfrac res;	//Contient l'ensemble des A_n-1 et l'ensemble Q_n
	mpz_inits(res.N, res.k, res.g, NULL);

	//Initialisation des variables
	mpz_t* A = (mpz_t*)malloc((rang+1)*sizeof(mpz_t));		//Le tableau contenant les A_n-1
	mpz_t* Q = (mpz_t*)malloc((rang+2)*sizeof(mpz_t));		//Le tableau contenant les Q_n
	mpz_t* P = (mpz_t*)malloc((rang+1)*sizeof(mpz_t));		//Éléments reliés au Q_n
	mpz_t* r = (mpz_t*)malloc((rang+1)*sizeof(mpz_t));		//Interviennent dans le calcul des A_n-1 & Q_n
	mpz_t* q = (mpz_t*)malloc(rang*sizeof(mpz_t));			//Idem
	mpz_t g, tempz;											//Idem, tempz = variable à tout faire	
	mpf_t sqrtkN, tempf, tempf2;									//sqrtkN = sqrt(k*N), tempf = variable à tout faire

	

	//Valeurs d'initialisation de la boucle
	mpz_inits(g, A[0], Q[0], r[0], tempz, NULL);
	mpf_inits(tempf, sqrtkN, tempf2, NULL);

	//Initialisation des différentes valeurs "indépendantes"

	mpz_set(tempz, N);
	mpz_mul(tempz, tempz, k);
	mpf_set_z(sqrtkN, tempz);
	mpf_sqrt(sqrtkN, sqrtkN);
	mpf_floor(tempf, sqrtkN);
	mpz_set_f(g, tempf);				//g = [sqrt(kN)]

	mpz_set(Q[0], k);
	mpz_mul(Q[0], Q[0], N);				//Q_-1 = kN = Q[0]
	mpz_set(r[0], g);					//r_-1 = r[0]
	mpz_set_ui(A[0], 1);				//A_-1 = A[0]

	
	//Calcul de P_0 & Q_0
	mpz_init_set_ui(Q[1], 1);			//Q_0 = Q[1]
	mpz_init_set_ui(P[0], 0);			//P_0 = P[0]	

	for(long long int i = 0; i < rang; i++)
	{
		switch(i)
		{
			case 0:

				//Calcul de q_0
				mpz_init_set(q[0], g);				//q_0 = [(sqrt(kN) + P_0)/Q_0] avec P_0 = 0 et Q_0 = 1

				//Calcul de A_0
				mpz_init_set(A[1],A[0]);
				mpz_mul(A[1], A[1], q[0]);			//A_0 = q_0*A_-1
				mpz_mod(A[1], A[1], N);				//On réduit mod N

				//Calcul de r_0
				mpz_init_set_ui(r[1], 0);			//r_0 = P_0 + g - q_0.Q_0 = 0 + g - g.1 = 0

				//Calcul de P_1
				mpz_init_set(P[1], g);				//P_1 = g - r_0 = g - 0 = g

				//Calcul de Q_1 = Q[2]
				mpz_init_set(Q[2], r[1]);		
				mpz_sub(Q[2], Q[2], r[0]);
				mpz_mul(Q[2], Q[2], q[0]);
				mpz_add(Q[2], Q[2], Q[0]);
				
				break;
			default:
				//Calcul q_i
				mpz_init(q[i]);	
				mpf_set_z(tempf, P[i]);
				mpf_set_z(tempf2, Q[i+1]);
				mpf_add(tempf, tempf, sqrtkN);				//sqrt(kN) + P_i
				mpf_div(tempf, tempf, tempf2);		
				mpf_floor(tempf, tempf);					//floor((sqrt(kN) + P_i)/Q_i)
				mpz_set_f(q[i], tempf);

				//Calcul de r_n = r[n+1]	
				mpz_init(r[i+1]);
				mpz_submul(r[i+1], q[i], Q[i+1]);
				mpz_add(r[i+1], r[i+1], P[i]);
				mpz_add(r[i+1], r[i+1], g);

				//Calcul de A_n = A[n+1] 		
				mpz_init_set(A[i+1],A[i]);
				mpz_mul(A[i+1], A[i+1], q[i]);		//A_i-1*q_i
				mpz_add(A[i+1], A[i+1], A[i-1]);	//A_i-1*q_i + A_i-2
				mpz_mod(A[i+1], A[i+1], N);			//réduction modulo N

				//Calcul P_n+1 
				mpz_init_set(P[i+1], g);
				mpz_sub(P[i+1], P[i+1], r[i+1]);		//P[n+1] = g - r_n = g - r[n+1]

				//Calcul Q_n+1 = Q[n+2]
				mpz_init_set(Q[i+2], r[i+1]); 
				mpz_sub(Q[i+2],	Q[i+2], r[i]);		//(r_n - r_n-1)
				mpz_mul(Q[i+2], Q[i+2], q[i]);		//q_n(r_n - r_n-1)
				mpz_add(Q[i+2], Q[i+2], Q[i]);		//Q_n-1 + q_n(r_n - r_n-1)
				
				break;
		}

	}

	//Test de routine pour voir si le développement de la fraction continue s'est bien passé
	mpz_t tempsqrt;
	mpz_init(tempsqrt);
	mpz_set(tempsqrt, k);
	mpz_mul(tempsqrt, tempsqrt, N);
	mpz_sqrt(tempsqrt, tempsqrt);
	mpz_mul_ui(tempsqrt, tempsqrt, 2);

	for(int i = 1; i < rang; i++)				//Éviter de commencer à i = 0 puisque cela représente Q_-1 qui n'intervient uniquement dans l'algo et non dans le développement en fraction continue
	{
		if(mpz_cmp(Q[i], tempsqrt) >= 0)	//Si Q_n >= 2sqrt(kN)
		{
			res.rang = 0;

			mpz_clears(g, tempz, NULL);
			mpf_clears(sqrtkN, tempf, tempf2, NULL);
			return res;
		}
	}

	//Assignation des tableaux dans le résultat
	res.A = A;
	res.Q = Q;
	res.r = r;
	res.q = q;
	res.P = P;
	mpz_set(res.N, N);
	mpz_set(res.k, k);
	res.rang = rang;
	mpz_set(res.g, g);

	//Libération de mémoire
	mpz_clears(g, tempz, NULL);
	mpf_clears(sqrtkN, tempf, tempf2, NULL);

	return res;
}
