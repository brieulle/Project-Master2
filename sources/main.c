#include"../entetes/lib.h"



int main(int argc, char** argv)//argv[1] = N, argv[2] = k, argv[3] = upperbound, argv[4] = rang, argv[5] = uplevel
{
	long long int upperbound = atoll(argv[3]);	//Borne supérieur de la base de factorisation
	long long int rang = atoll(argv[4]);		//Rang jusqu'auquel on développera la fraction continue
	long long int uplevel = atoll(argv[5]);	//Entier par lequel on augementera la taille de du rang et de la borne pour la base de factorisation

	mpz_t N, k;								//N entier à factoriser, k entier servant à "élargir" le développement de N en fraction continue
	mpz_inits(N, k, NULL);

	mpz_set_str(N,argv[1],10);
	mpz_set_ui(k,atoi(argv[2]));

	boucleFactorisation(N, k, upperbound, rang, uplevel);

	mpz_clears(N, k, NULL);
	return 0;
}
