#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define CONST_Incidenza			50
#define CONST_Seguiti			50
#define CONST_Seguaci			1
#define CONST_SeguitiPesati		100
#define CONST_Sorensen			5
#define CONST_HDI			1
#define CONST_RA			100
#define CONST_RWR			100
#define CONST_RWRP			100

static const double Crwr = 0.5;
static const int numIterazioniPseudoInversa = 3;

static const int hashA = 1;
static const int hashB = 1;
static const int prime = 1000000007;
static const int DIMENSIONETABELLA = 1000;

//INFO HASH
// con tabella hash grande 1000 si hanno sovrapposizioni di al più 5 elementi in uno stesso bucket
//INFO GRAFO
// sono 959 nodi, 5000 archi
// max archi in uscita = 94
// max archi in entrata = 22

//TIPI
// tipo 0: rank 1 - 102, 			sono 43 elementi
// tipo 1: rank 104 - 7196, 		sono 127 elementi
// tipo 2: rank 7236 - 96059		sono 255 elementi
// tipo 3: rank 96228 - 2147483647 	sono 429 elementi
// senza tipo						sono 105 elementi


typedef struct nodo {
	char *url;
	int type;
	long long int rank;
	struct lista *inizioListaAdiacentiSeguaci;		//siti che hanno copiato il nodo
	struct lista *inizioListaAdiacentiSeguiti;		//siti da cui il nodo ha copiato
	int numSeguaci;
	int numSeguiti;
	int indice;										//indice nell'array dei nodi
} nodo_t;
typedef struct arco {
	int nodo1;
	int nodo2;
	double valore;
} arco_t;
typedef struct lista {
	struct lista *next;
	nodo_t *el;
} lista_t;
typedef struct hashTable {
    int a, b;
    int n_bucket;
	int *overlap;
    lista_t **table;
} hashTable_t;


nodo_t* newNodo (char *url) {
	nodo_t *new = malloc(sizeof(nodo_t));
	new->url = malloc((1+strlen(url))*sizeof(char));
	stpcpy(new->url,url);
	new->type = 2;
	new->rank = 0;
	new->inizioListaAdiacentiSeguaci = NULL;
	new->inizioListaAdiacentiSeguiti = NULL;
	new->numSeguaci = 0;
	new->numSeguiti = 0;
	return new;
}
lista_t* newElementoLista (nodo_t* elemento, lista_t* next) {
	lista_t* new = malloc(sizeof(lista_t*));
	new->el = elemento;
	new->next = next;
	return new;
}
int comparaArchi (const void * elem1, const void * elem2) {
    arco_t* x1 = (arco_t*) elem1;
    arco_t* x2 = (arco_t*) elem2;
    if (x1->valore > x2->valore) return -1;
    if (x1->valore < x2->valore) return 1;
    return 0;
}
void stampaNodo(nodo_t *sito) {
	lista_t* provv;
	nodo_t* nodo;
	printf("Il sito %s  con rank %lli  ", sito->url, sito->rank);
	if (sito->type!=-1) printf("di tipo %d\n", sito->type);
	else printf ("senza tipo\n");
	printf("\t ha %d followers\n", sito->numSeguaci);
		for (provv=sito->inizioListaAdiacentiSeguaci; provv!=NULL; provv=provv->next) {
			nodo = provv->el;
			printf("\t\t%s(%lli)\n",nodo->url,nodo->rank);
		}
	printf("\t ha %d copiati\n", sito->numSeguiti);
		for (provv=sito->inizioListaAdiacentiSeguiti; provv!=NULL; provv=provv->next) {
			nodo = provv->el;
			printf("\t\t%s(%lli)\n",nodo->url,nodo->rank);
		}
}


//OPERAZIONI TABELLA HASH
hashTable_t* newHashTable (int expected_size) {
	hashTable_t* new = malloc(sizeof(hashTable_t));
	new->a = hashA;
	new->b = hashB;
	new->n_bucket = expected_size;
	new->overlap = calloc(expected_size, sizeof(int));
	new->table = calloc(expected_size, sizeof(lista_t*));
    return new;
}
int hash (hashTable_t* ht, char* parola) {
	unsigned long long numero = 0;
	int i, lunghezzaStringa = strlen(parola);
	for (i = 0; i < lunghezzaStringa-2; i++) numero = 27*(numero + 1 + parola[i] - 'a');
    return ((((  (ht->a)*numero +  (ht->b)  ) % prime) + prime ) % prime) % (ht->n_bucket);
}
nodo_t* insert (hashTable_t* ht, char* url) {
	int h = hash(ht,url);
	(ht->overlap[h])++;
	lista_t *provv = (ht->table)[h];
	
	if (provv == NULL) {
		ht->table[h] = newElementoLista(newNodo(url), NULL);
		return ht->table[h]->el;
	}
	else {
		while (provv->next!=NULL) 	provv=provv->next;
		provv->next = newElementoLista(newNodo(url), NULL);
		return provv->next->el;
	}
}
nodo_t* urlToPointer (hashTable_t* ht, char *url) {
	lista_t *provv;
	nodo_t *nodo;
	for (provv=(ht->table)[hash(ht,url)]; provv!=NULL; provv=provv->next) {
		nodo = provv->el;
		if (strcmp(nodo->url, url)==0) return nodo;
	}
	return NULL;
}


//OPERAZIONI CON MATRICI
void multiplicationMatrix(double* a, double* b, double* c, int size) {
	int i,j,k;
	for (i=0;i<size;i++) for (j=0;j<size;j++) {
		c[i*size+j] = 0;
		for (k=0;k<size;k++) c[i*size+j] += a[i*size+k] * b[k*size+j];
	}
}
void multiplicationMatrixPerScalar(double* a, double scalar, int size) {
	int i,j;
	for (i=0;i<size;i++) for (j=0;j<size;j++) 	a[i*size+j] *= scalar;
}
void transposeMatrix(double* a, int size) {
	double temp;
	int i,j;
	for (i=0;i<size;i++) for (j=0;j<i;j++) {
		temp = a[i*size+j];
		a[i*size+j] = a[j*size+i];
		a[j*size+i] = temp;
	}
}
double* pseudoInverseMatrix(double* P, int numIterazioni, int size) {
	//calcola l'inversa di Id - P tramite serie geometrica troncata
	double* potenzaDiP = malloc(size*size*sizeof(double));
	double* potenzaDiPsuccessiva = malloc(size*size*sizeof(double));
	int i,j,k;	
	double* result = calloc(size*size,sizeof(double));
	for (i=0;i<size;i++) result[i*size+i] = 1;
	for (i=0;i<size;i++) for (j=0;j<size;j++)	result[i*size+j] += (potenzaDiP[i*size+j] = P[i*size+j]);
	for (k=1;k<numIterazioni;k++) {
		multiplicationMatrix(P, potenzaDiP, potenzaDiPsuccessiva, size);
		for (i=0;i<size;i++) for (j=0;j<size;j++)	result[i*size+j] += (potenzaDiP[i*size+j] = potenzaDiPsuccessiva[i*size+j]);
	}
	return result;
}
void normalizza(double* a, int size) {
	int i,j;
	double max=0;
	for(i=0;i<size;i++) for(j=0;j<size;j++)	if (a[i*size+j] > max)	max=a[i*size+j];
	if (max!=0)	for(i=0;i<size;i++) for(j=0;j<size;j++)	a[i*size+j] /= max;
}


//CHECK "TRANSITIVITA" DEI RISULTATI
double transitivityCheck_CasoMediano(double* a, double range, int size) {
	// A[i,j] deve stare nel range giusto
	int i,j,k,coppieBuone=0,coppieCattive=0;
	for(i=0;i<size;i++)		for(j=0;j<i;j++)	for(k=0;k<size;k++)		if(k!=i && k!=j) {
			if ((a[i*size+j] > a[i*size+k]*a[k*size+j]-range) && (a[i*size+j] < a[i*size+k]*a[k*size+j]+range))		coppieBuone++;
			else																									coppieCattive++;
	}
	return ((double)coppieBuone)/(coppieBuone+coppieCattive);
}
double transitivityCheck_CasoMinorato(double* a, int size) {
	// A[i,j] deve valere almeno tot
	int i,j,k,coppieBuone=0,coppieCattive=0;
	for(i=0;i<size;i++)		for(j=0;j<i;j++)	for(k=0;k<size;k++)		if(k!=i && k!=j) {
			if ( a[i*size+j]  >=  1 - (1-a[i*size+k]) - (1-a[k*size+j]) )	coppieBuone++;
			else													coppieCattive++;
	}
	return ((double)coppieBuone)/(coppieBuone+coppieCattive);
}
double transitivityCheck_CasoMaggiorato(double* a, int size) {
	// A[i,j] deve valere al massimo tot
	int i,j,k,differenza,coppieBuone=0,coppieCattive=0;
	for(i=0;i<size;i++)		for(j=0;j<i;j++)	for(k=0;k<size;k++)		if(k!=i && k!=j) {
			differenza = a[i*size+k] - a[k*size+j];
			if (differenza < 0) differenza = -differenza;
			if ( a[i*size+j] <= 1 - differenza )	coppieBuone++;
			else									coppieCattive++;
	}
	return ((double)coppieBuone)/(coppieBuone+coppieCattive);
}

int main() {
	FILE *archi = fopen("Input - archi.txt","r");
	FILE *ranks = fopen("Input - ranks.txt","r");
	int numNodi=0;
	int i, j, c1, type;
	long long int rank;
	float c2;
	double media,numero,max;
	char urlStart[50];
	char urlEnd[50];
	nodo_t *nodoStart, *nodoEnd, *nodo1, *nodo2, *nodo;
	lista_t *provv, *provv1, *provv2;
	hashTable_t* ht = newHashTable(DIMENSIONETABELLA);
	
	
	//METTO NODI E ARCHI IN TABELLA HASH
	for (i=0; i<5000; i++) {
		fscanf(archi, "%d %s %s %d %f %f %f", &c1, urlStart, urlEnd, &c1, &c2, &c2, &c2);
		//INSERISCO I NODI SE NON CI SONO GIA
		if ( (nodoStart = urlToPointer(ht, urlStart)) == NULL)	{			
			numNodi++;
			nodoStart = insert(ht,urlStart);
		}
		if ( (nodoEnd = urlToPointer(ht, urlEnd)) == NULL)	{
			numNodi++;
			nodoEnd = insert(ht,urlEnd);
        }
		//INSERISCO L'ARCO	
		nodoStart->inizioListaAdiacentiSeguiti = newElementoLista(nodoEnd, nodoStart->inizioListaAdiacentiSeguiti);
		nodoStart->numSeguiti++;
		nodoEnd->inizioListaAdiacentiSeguaci = newElementoLista(nodoStart, nodoEnd->inizioListaAdiacentiSeguaci);
		nodoEnd->numSeguaci++;
    }
	printf("\nCi sono %d nodi, con tabella hash grande %d\n",numNodi,DIMENSIONETABELLA);

	//AGGIUNGO TIPO E RANK DATI DA AMAZON ALEXA
	for (i=0; i<854; i++) {
		fscanf(ranks, "%d %s %d %lli",&c1, urlStart, &type, &rank);
		nodo = urlToPointer(ht, urlStart);
		nodo->type = type;
		nodo->rank = rank;
    }

	//METTO I NODI IN UN ARRAY (di puntatori)
	j=0;
	nodo_t** arrayNodi = malloc(numNodi*sizeof(nodo_t*));
	for (i=0; i < ht->n_bucket; i++)
		for (provv = ht->table[i]; provv!=NULL; provv=provv->next) {
			nodo = provv->el;
			nodo->indice = j;
			arrayNodi[j] = nodo;
			j++;
	}

	
	
	
	
	
	
	//CALCOLO ALCUNI INDICI DI (supposta) SIMILITUDINE
	
	//INDICE 0
	//TABELLA INCIDENZA
	int *tabellaIncidenza = calloc(numNodi*numNodi, sizeof(int));
	double *tabellaIncidenzaSimmetrica = malloc(numNodi*numNodi*sizeof(double));
	if (CONST_Incidenza!=0) {
		for (i=0;i<numNodi;i++)
			for (provv = arrayNodi[i]->inizioListaAdiacentiSeguiti; provv!=NULL; provv=provv->next) {
				nodo = provv->el;
				tabellaIncidenza[i*numNodi + nodo->indice] = 1;
		}
		printf("\nHo calcolato la tabella di incidenza\n");
		//simmetrizzo
		for (i=0;i<numNodi;i++) for (j=0;j<i;j++) tabellaIncidenzaSimmetrica[i*numNodi+j] = tabellaIncidenza[i*numNodi+j]+tabellaIncidenza[j*numNodi+i];
	}
	
	
	//INDICE 1 (CN)
	//NUMERO COPIATI IN COMUNE
	double *seguitiInComune = calloc(numNodi*numNodi, sizeof(double));
	FILE *seguitiInComuneFILE;
	if (CONST_Seguiti!=0 && (CONST_Sorensen!=0 || CONST_HDI!=0)) {
		for (i=0;i<numNodi;i++)
			for (provv1=arrayNodi[i]->inizioListaAdiacentiSeguaci; provv1!=NULL; provv1=provv1->next) {
				nodo1 = provv1->el;
				for (provv2=provv1->next; provv2!=NULL; provv2=provv2->next) {
					nodo2 = provv2->el;
					if (nodo1->indice > nodo2->indice)	seguitiInComune[nodo1->indice * numNodi + nodo2->indice]++;
					else								seguitiInComune[nodo2->indice * numNodi + nodo1->indice]++;
				}
		}
		seguitiInComuneFILE = fopen("Index 1 - Common Neighbors Seguiti.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(seguitiInComuneFILE,"%0.0f ", seguitiInComune[i*numNodi+j]);
			fprintf(seguitiInComuneFILE,"\n");
		}
		printf("Ho calcolato Common Neighbour Index (seguiti)\n");
	}
	//NUMERO SEGUACI IN COMUNE
	FILE *seguaciInComuneFILE;
	double *seguaciInComune = calloc(numNodi*numNodi, sizeof(double));
	if (CONST_Seguaci!=0) {
		for (i=0;i<numNodi;i++)
			for (provv1=arrayNodi[i]->inizioListaAdiacentiSeguiti; provv1!=NULL; provv1=provv1->next) {
				nodo1 = provv1->el;
				for (provv2=provv1->next; provv2!=NULL; provv2=provv2->next) {
					nodo2 = provv2->el;
					if (nodo1->indice > nodo2->indice)	seguaciInComune[nodo1->indice * numNodi + nodo2->indice]++;
					else								seguaciInComune[nodo2->indice * numNodi + nodo1->indice]++;
				}
		}
		seguaciInComuneFILE = fopen("Index 1 - Common Neighbors Seguaci.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(seguaciInComuneFILE,"%0.0f ", seguaciInComune[i*numNodi+j]);
			fprintf(seguaciInComuneFILE,"\n");
		}
		printf("Ho calcolato Common Neighbour Index (seguaci)\n");
	}
	
	
	//INDICE 2
	//NUMERO COPIATI IN COMUNE PESATI CON RANK		(Se due siti copiano entrambi un sito semisconosciuto vuol dire che avevano un buon motivo per farlo)
	double *seguitiInComunePESATI = calloc(numNodi*numNodi, sizeof(double));
	FILE *seguitiInComunePesatiFILE;
	if (CONST_SeguitiPesati!=0) {
		for (i=0;i<numNodi;i++)
			for (provv1=arrayNodi[i]->inizioListaAdiacentiSeguaci; provv1!=NULL; provv1=provv1->next) {
				nodo1 = provv1->el;
				for (provv2=provv1->next; provv2!=NULL; provv2=provv2->next) {
					nodo2 = provv2->el;
					if (nodo1->indice > nodo2->indice)	seguitiInComunePESATI[nodo1->indice * numNodi + nodo2->indice] += arrayNodi[i]->type;
					else 								seguitiInComunePESATI[nodo2->indice * numNodi + nodo1->indice] += arrayNodi[i]->type;
				}
		}
		seguitiInComunePesatiFILE = fopen("Index 2 - Common Neighbors Seguiti Pesati.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(seguitiInComunePesatiFILE,"%0.0f ", seguitiInComunePESATI[i*numNodi+j]);
			fprintf(seguitiInComunePesatiFILE,"\n");
		}
		printf("Ho calcolato Common Neighbour Pesati Index\n");
	}
	
	
	//INDICE 3		(Sorensen index)
	//NUMERO COPIATI IN COMUNE DIVISO PER MEDIA DI NUMERO DI NODI COPIATI DAI DUE
	double *seguitiInComuneSorensen = calloc(numNodi*numNodi, sizeof(double));
	FILE *seguitiInComuneSorensenFILE;
	if (CONST_Sorensen!=0) {
		for (i=0;i<numNodi;i++) for(j=0;j<i;j++) {
			media = ( arrayNodi[i]->numSeguiti + arrayNodi[j]->numSeguiti )/2.0;
			if (media!=0)	seguitiInComuneSorensen[i*numNodi+j] = seguitiInComune[i*numNodi+j] / media;
			else 			seguitiInComuneSorensen[i*numNodi+j] = 0;
		}
		seguitiInComuneSorensenFILE = fopen("Index 3 - Sorensen.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(seguitiInComuneSorensenFILE,"%0.3f ", seguitiInComuneSorensen[i*numNodi+j]);
			fprintf(seguitiInComuneSorensenFILE,"\n");
		}
		printf("Ho calcolato Sorensen Index\n");
	}

	
	//INDICE 4		(HDI)
	//NUMERO COPIATI IN COMUNE DIVSO PER MAX DI NUMERO DI NODI COPIATI DAI DUE
	double *seguitiInComuneHubDepressed = calloc(numNodi*numNodi, sizeof(double));
	FILE *seguitiInComuneHubDepressedFILE;
		if (CONST_HDI!=0) {
		for (i=0;i<numNodi;i++) for(j=0;j<i;j++) {
			max = arrayNodi[i]->numSeguiti;
			if (arrayNodi[j]->numSeguiti > max)		max = arrayNodi[j]->numSeguiti;
			if (max!=0)		seguitiInComuneHubDepressed[i*numNodi+j] = seguitiInComune[i*numNodi+j] / max;
			else 			seguitiInComuneHubDepressed[i*numNodi+j] = 0;
		}
		seguitiInComuneHubDepressedFILE = fopen("Index 4 - HDI.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(seguitiInComuneHubDepressedFILE,"%0.3f ", seguitiInComuneHubDepressed[i*numNodi+j]);
			fprintf(seguitiInComuneHubDepressedFILE,"\n");
		}
		printf("Ho calcolato HDI Index\n");
	}
	
	
	//INDICE 5		(RA)
	//SOMMA SUI SEGUITI IN COMUNE DI ( 1 / NUMERO SEGUACI )
	double *resourceAllocation = calloc(numNodi*numNodi, sizeof(double));
	FILE *resourceAllocationFILE;
	if (CONST_RA!=0) {
		for (i=0;i<numNodi;i++)		//per comodità inverto l'ordine delle sommatorie e sommo prima sui nodi
			for (provv1=arrayNodi[i]->inizioListaAdiacentiSeguaci; provv1!=NULL; provv1=provv1->next) {
				nodo1 = provv1->el;
				for (provv2=provv1->next; provv2!=NULL; provv2=provv2->next) {
					nodo2 = provv2->el;
					if (nodo1->indice > nodo2->indice)	resourceAllocation[nodo1->indice * numNodi + nodo2->indice] += 1.0/arrayNodi[i]->numSeguaci;
					else 								resourceAllocation[nodo2->indice * numNodi + nodo1->indice] += 1.0/arrayNodi[i]->numSeguaci;
				}
		}
		resourceAllocationFILE = fopen("Index 5 - RA.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(resourceAllocationFILE,"%0.3f ", resourceAllocation[i*numNodi+j]);
			fprintf(resourceAllocationFILE,"\n");
		}
		printf("Ho calcolato RA Index\n");
	}
	
	
	//INDICE 6		(RWR)
	//RANDOM WALK WITH RESTART
	double* randomWalkWithRestart = calloc(numNodi*numNodi,sizeof(double));
	double* transitionMatrix;
	double* pseudoInversa;
	double* identity;
	FILE *randomWalkWithRestartFILE;
	FILE *identityFILE;
	if (CONST_RWR!=0) {
		//calcolo la matrice di transizione
		transitionMatrix = calloc(numNodi*numNodi,sizeof(double));
		for (i=0;i<numNodi;i++) {
			if (arrayNodi[i]->numSeguiti != 0) {
				numero = 1.0 / arrayNodi[i]->numSeguiti;
				for (provv=arrayNodi[i]->inizioListaAdiacentiSeguiti; provv!=NULL; provv=provv->next)	transitionMatrix[i*numNodi + provv->el->indice] = numero;
			}
			else {
				numero = 1.0 / (numNodi-1);
				for (j=0;j<numNodi;j++)		if (i!=j)	transitionMatrix[i*numNodi+j] = numero;
			}
		}
		//calcolo la pseudoInversa = (Id - c*transitionMatrixT)^(-1)
		transposeMatrix(transitionMatrix, numNodi);
		multiplicationMatrixPerScalar(transitionMatrix, Crwr, numNodi);
		pseudoInversa = pseudoInverseMatrix(transitionMatrix, numIterazioniPseudoInversa, numNodi);
		//check inversa
		for (i=0;i<numNodi;i++) for(j=0;j<numNodi;j++) transitionMatrix[i*numNodi+j] = -transitionMatrix[i*numNodi+j];
		for (i=0;i<numNodi;i++) transitionMatrix[i*numNodi+i] += 1;
		identity = calloc(numNodi*numNodi,sizeof(double));
		multiplicationMatrix(transitionMatrix, pseudoInversa, identity, numNodi);
		identityFILE = fopen("Identity Check.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<numNodi;j++) fprintf(identityFILE,"%5.5f ", identity[i*numNodi+j]);
			fprintf(identityFILE,"\n");
		}
		//calcolo il Random Walk With Restart Index
		multiplicationMatrixPerScalar(pseudoInversa, 1-Crwr, numNodi);
		for (i=0;i<numNodi;i++) 	for (j=0;j<i;j++)	randomWalkWithRestart[i*numNodi+j] = pseudoInversa[i*numNodi+j] + pseudoInversa[j*numNodi+i];
		randomWalkWithRestartFILE = fopen("Index 6 - RWR.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(randomWalkWithRestartFILE,"%0.3f ", randomWalkWithRestart[i*numNodi+j]);
			fprintf(randomWalkWithRestartFILE,"\n");
		}
		printf("Ho calcolato il RWR Index\n");
	}
	
	
	//INDICE 7		(RWR pesato)
	//RANDOM WALK WITH RESTART CON MATRICE DI TRANSIZIONE PESATA CON RANK
	int sommaTypeSeguiti;
	double* randomWalkWithRestartPesato = calloc(numNodi*numNodi,sizeof(double));
	double* transitionMatrixPesata;
	double* pseudoInversaPesata;
	double* identityPesata;
	FILE *randomWalkWithRestartPesatoFILE;
	FILE *identityPesataFILE;
	if (CONST_RWRP!=0) {
		//calcolo la matrice di transizione
		transitionMatrixPesata = calloc(numNodi*numNodi,sizeof(double));
		for (i=0;i<numNodi;i++) {
			if (arrayNodi[i]->numSeguiti != 0) {
				sommaTypeSeguiti=0;
				for (provv=arrayNodi[i]->inizioListaAdiacentiSeguiti; provv!=NULL; provv=provv->next)	sommaTypeSeguiti += provv->el->type + 1;
				for (provv=arrayNodi[i]->inizioListaAdiacentiSeguiti; provv!=NULL; provv=provv->next) {
					nodo = provv->el;
					transitionMatrixPesata[i*numNodi + nodo->indice] = (nodo->type + 1.0)/sommaTypeSeguiti;
				}
			}
			else {
				numero = 1.0 / (numNodi-1);
				for (j=0;j<numNodi;j++)		if (i!=j)	transitionMatrix[i*numNodi+j] = numero;
			}
		}
		//calcolo la pseudoInversa = (Id - c*transitionMatrixT)^(-1)
		transposeMatrix(transitionMatrixPesata, numNodi);
		multiplicationMatrixPerScalar(transitionMatrixPesata, Crwr, numNodi);
		pseudoInversaPesata = pseudoInverseMatrix(transitionMatrixPesata, numIterazioniPseudoInversa, numNodi);
		//check inversa
		for (i=0;i<numNodi;i++) for(j=0;j<numNodi;j++) transitionMatrixPesata[i*numNodi+j] = -transitionMatrixPesata[i*numNodi+j];
		for (i=0;i<numNodi;i++) transitionMatrixPesata[i*numNodi+i] += 1;
		identityPesata = calloc(numNodi*numNodi,sizeof(double));
		multiplicationMatrix(transitionMatrixPesata, pseudoInversaPesata, identityPesata, numNodi);
		identityPesataFILE = fopen("Identity Check PESATA.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<numNodi;j++) fprintf(identityPesataFILE,"%5.5f ", identityPesata[i*numNodi+j]);
			fprintf(identityPesataFILE,"\n");
		}
		//calcolo il Random Walk With Restart Index
		multiplicationMatrixPerScalar(pseudoInversaPesata, 1-Crwr, numNodi);
		for (i=0;i<numNodi;i++) 	for (j=0;j<i;j++)	randomWalkWithRestartPesato[i*numNodi+j] = pseudoInversaPesata[i*numNodi+j] + pseudoInversaPesata[j*numNodi+i];
		randomWalkWithRestartPesatoFILE = fopen("Index 7 - RWRP.txt","w");
		for (i=0;i<numNodi;i++) {
			for (j=0;j<i;j++) fprintf(randomWalkWithRestartPesatoFILE,"%0.3f ", randomWalkWithRestartPesato[i*numNodi+j]);
			fprintf(randomWalkWithRestartPesatoFILE,"\n");
		}
		printf("Ho calcolato il RWRP Index\n");
	}
	
	
	
	//INDICE FINALE
	//Somma di tutti gli indici pesati con la costante scelta
	normalizza(tabellaIncidenzaSimmetrica, numNodi);
	normalizza(seguitiInComune, numNodi);
	normalizza(seguaciInComune, numNodi);
	normalizza(seguitiInComunePESATI, numNodi);
	normalizza(seguitiInComuneSorensen, numNodi);
	normalizza(seguitiInComuneHubDepressed, numNodi);
	normalizza(resourceAllocation, numNodi);
	normalizza(randomWalkWithRestart, numNodi);
	normalizza(randomWalkWithRestartPesato, numNodi);
	
	printf("\n\nIndice \tUpperBound \tLowerBound \tCheck 0.25 \tCheck 0.1 \tCheck 0.01\n");
	printf ("I. Simmetrica =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n",
													transitivityCheck_CasoMaggiorato(tabellaIncidenzaSimmetrica,numNodi),
													transitivityCheck_CasoMinorato(tabellaIncidenzaSimmetrica,numNodi),
													transitivityCheck_CasoMediano(tabellaIncidenzaSimmetrica,0.25,numNodi),
													transitivityCheck_CasoMediano(tabellaIncidenzaSimmetrica,0.1,numNodi),
													transitivityCheck_CasoMediano(tabellaIncidenzaSimmetrica,0.01,numNodi) );
	printf ("SeguitiIC =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 		
													transitivityCheck_CasoMaggiorato(seguitiInComune,numNodi),
													transitivityCheck_CasoMinorato(seguitiInComune,numNodi),
													transitivityCheck_CasoMediano(seguitiInComune,0.25,numNodi),
													transitivityCheck_CasoMediano(seguitiInComune,0.1,numNodi),
													transitivityCheck_CasoMediano(seguitiInComune,0.01,numNodi) );
	printf ("SeguaciIC =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 		
													transitivityCheck_CasoMaggiorato(seguaciInComune,numNodi),
													transitivityCheck_CasoMinorato(seguaciInComune,numNodi),
													transitivityCheck_CasoMediano(seguaciInComune,0.25,numNodi),
													transitivityCheck_CasoMediano(seguaciInComune,0.1,numNodi),
													transitivityCheck_CasoMediano(seguaciInComune,0.01,numNodi) );
	printf ("SeguitiICP =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 		 
													transitivityCheck_CasoMaggiorato(seguitiInComunePESATI,numNodi),
													transitivityCheck_CasoMinorato(seguitiInComunePESATI,numNodi),
													transitivityCheck_CasoMediano(seguitiInComunePESATI,0.25,numNodi),
													transitivityCheck_CasoMediano(seguitiInComunePESATI,0.1,numNodi),
													transitivityCheck_CasoMediano(seguitiInComunePESATI,0.01,numNodi) );
	printf ("SegSorensen =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 	
													transitivityCheck_CasoMaggiorato(seguitiInComuneSorensen,numNodi),
													transitivityCheck_CasoMinorato(seguitiInComuneSorensen,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneSorensen,0.25,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneSorensen,0.1,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneSorensen,0.01,numNodi) );
	printf ("HDI =\t\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 			 
													transitivityCheck_CasoMaggiorato(seguitiInComuneHubDepressed,numNodi),
													transitivityCheck_CasoMinorato(seguitiInComuneHubDepressed,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneHubDepressed,0.25,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneHubDepressed,0.1,numNodi),
													transitivityCheck_CasoMediano(seguitiInComuneHubDepressed,0.01,numNodi) );
	printf ("RA =\t\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 			 
													transitivityCheck_CasoMaggiorato(resourceAllocation,numNodi),
													transitivityCheck_CasoMinorato(resourceAllocation,numNodi),
													transitivityCheck_CasoMediano(resourceAllocation,0.25,numNodi),
													transitivityCheck_CasoMediano(resourceAllocation,0.1,numNodi),
													transitivityCheck_CasoMediano(resourceAllocation,0.01,numNodi) );
	printf ("RWR =\t\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 			
													transitivityCheck_CasoMaggiorato(randomWalkWithRestart,numNodi),
													transitivityCheck_CasoMinorato(randomWalkWithRestart,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestart,0.25,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestart,0.1,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestart,0.01,numNodi) );
	printf ("RWRP =\t\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 			
													transitivityCheck_CasoMaggiorato(randomWalkWithRestartPesato,numNodi),
													transitivityCheck_CasoMinorato(randomWalkWithRestartPesato,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestartPesato,0.25,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestartPesato,0.1,numNodi),
													transitivityCheck_CasoMediano(randomWalkWithRestartPesato,0.01,numNodi) );
	
	double *similarityIndex = calloc(numNodi*numNodi, sizeof(double));
	for (i=0;i<numNodi;i++) for(j=0;j<i;j++) {
		similarityIndex[i*numNodi+j] = 
				CONST_Incidenza * tabellaIncidenzaSimmetrica[i*numNodi+j] +
				CONST_Seguiti * seguitiInComune[i*numNodi+j] +
				CONST_Seguaci * seguaciInComune[i*numNodi+j] +
				CONST_SeguitiPesati * seguitiInComunePESATI[i*numNodi+j] +
				CONST_Sorensen * seguitiInComuneSorensen[i*numNodi+j] +
				CONST_HDI * seguitiInComuneHubDepressed[i*numNodi+j] +
				CONST_RA * resourceAllocation[i*numNodi+j] +
				CONST_RWR * randomWalkWithRestart[i*numNodi+j] +
				CONST_RWRP * randomWalkWithRestartPesato[i*numNodi+j];
	}
	
	
	normalizza(similarityIndex, numNodi);
	printf ("\n\nFINALE =\t %0.3f,\t %0.3f,\t\t %0.5f,\t %0.5f,\t %0.5f\n", 		
													transitivityCheck_CasoMaggiorato(similarityIndex,numNodi),
													transitivityCheck_CasoMinorato(similarityIndex,numNodi),
													transitivityCheck_CasoMediano(similarityIndex,0.25,numNodi),
													transitivityCheck_CasoMediano(similarityIndex,0.1,numNodi),
													transitivityCheck_CasoMediano(similarityIndex,0.01,numNodi) );

	

	//ordino gli archi per valore finale
	arco_t* arrayArchi = malloc(numNodi*(numNodi-1)/2 * sizeof(arco_t));
	int counter=0;
	for (i=0;i<numNodi;i++) for (j=0;j<i;j++) {
		arrayArchi[counter].nodo1 = arrayNodi[i]->indice;
		arrayArchi[counter].nodo2 = arrayNodi[j]->indice;
		arrayArchi[counter].valore = similarityIndex[i*numNodi+j];
		counter++;
	}
	qsort(arrayArchi, counter, sizeof(arco_t), comparaArchi);
	FILE *OrderedPairsFILE = fopen("ORDERED PAIRS.txt","w");
	fprintf(OrderedPairsFILE, "Punteggio \tUrl 1 \t- Url 2\n");
	for (i=0; i<counter; i++) {
		fprintf(OrderedPairsFILE, "%5.6f \t%s \t- %s\n", arrayArchi[i].valore, arrayNodi[arrayArchi[i].nodo1]->url, arrayNodi[arrayArchi[i].nodo2]->url);
	}
	
	FILE *BindingIndiceUrlFILE = fopen("Binding Indice-Url.txt","w");
	fprintf(BindingIndiceUrlFILE, "Indice \tUrl\n");
	for (i=0; i<numNodi; i++) {
		fprintf(BindingIndiceUrlFILE, "%d\t- %s\n", i, arrayNodi[i]->url);
	}
	
	return 0;
}
