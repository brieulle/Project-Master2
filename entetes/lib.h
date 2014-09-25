#ifndef LIB_H
#define LIB_H

//entêtes de base
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<gmp.h>

/*Fonctions diverses*/

//La borne-ème fraction convergente (modifiée)
typedef struct
{
	mpz_t* A;
	mpz_t* Q;
	mpz_t* r;
	mpz_t* q;
	mpz_t* P;

	mpz_t N;		//L'entier à factoriser
	mpz_t k;		//Le coeff multiplicateur (au cas où la période de racine de N est trop courte
	long long unsigned int borne;
	mpz_t g;		//La partie entière de la racine carrée de k*N
} cfrac;

typedef struct
{
	mpz_t* tab;
	long long int taille;
//	unsigned int UpperBound;
} base; 

void cfrac_display(const cfrac);
void base_display(const base);

base gen_primeBase(const long long unsigned int, const mpz_t, const mpz_t, FILE*);	//UpperBound : borne jusqu'où on prend les nombres premiers, N, k; Fonction qui sert à générer la base de nombres premiers en fonction de N et k 

/*Fonctions de l'étape A*/

cfrac expand(const mpz_t, const mpz_t, const mpz_t); //N: l'entier dont la racine carrée est à développer en fraction continue, n_0 : le rang jusqu'où on développe la fraction, k le coefficient pour la racine carrée

/*Fonctions de l'étape B*/

typedef struct
{
	long long int* n; 			//Contient les indices des Q_n B-friable
	mpz_t* valQ;		//Tableau qui contiendra les valuations p-adique modulo 2 des Q_n pour p dans la base de factorisation, le vecteur est de la forme e_n = (a_0,...,a_r), où a_0 est le bit de parité de l'indice de n et a_i la puissance modulo deux 
	mpz_t* h;			//Le tableau qui contient les vecteurs historiques	
	long long unsigned int nb_Q;	//Le nombre de Q_n B-friable
	mpz_t* nQ;			//Les Q_n B-friable

} preS_set;

typedef struct
{
	mpz_t* Qs;				//Tableau contenant les Q^2 = produit signé de Q_n
	mpz_t* As;				//Tableau contenant
	long long unsigned int taille;	//Le nombre de pair (sert éventuellement de critère de réussite pour l'étape B
} S_set;


preS_set pairingSquare(mpz_t*, const long long unsigned int, const base, const mpz_t); //Q, borne, base de factorisation, l'entier N à factoriser; Fonction qui servira à calculer les couples A-S tel que décrit dans l'algorithme

preS_set factQ(mpz_t*, const long long unsigned int, const base, const mpz_t); //Q, borne, primebase, N; Fonction dont le but est de factoriser les Q_n selon la base de factorisation

/*Fonctions de l'étape C*/

S_set S_set_affect(const preS_set, mpz_t*, const mpz_t);	//F, le preset avec les Q_n factorisés, les A_n-1

void isFactored(const mpz_t, S_set, mpz_t);	//N l'entier à factoriser, S_set l'ensemble des pairs A - Q, le facteur (égal à -1 si échec)

/*Fonctions de Factorisation*/
int factorisation(const mpz_t, const mpz_t, const long long unsigned int, long unsigned int, FILE*); //N, k, borne, UpperBound, Fichier contenant les nombres premiers
void boucleFactorisation(const mpz_t, const mpz_t);
#endif
