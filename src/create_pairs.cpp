//
// Created by Adrián on 22/4/24.
//
#include "iostream"

extern "C" {
#include "k2-tree/matrix.h"
}

typedef struct
{ uint32_t s,p,o;
} triple;

#define step 1000000

int main(int argc, char **argv) {

    if (argc < 4) {
        std::cerr << "\tUsage: " << argv[0] << " <dataset> <n_preds> <n_triples>" << std::endl;
        exit(1);
    }

    std::string dataset = argv[1];
    uint n_preds = atoi(argv[2]);
    uint n_triples = atoi(argv[3]);

    uint S = n_preds+1;
    uint maxV = 0;
     triple *Set;
        uint j,i,n;
        FILE *f;
        FILE *files[1001];

        Set = (triple*) myalloc(n_triples*sizeof(triple));

        f = fopen(dataset.c_str(),"r");
        for (n=0;n<n_triples;n++)
        { fscanf(f,"%i %i %i\n",&Set[n].s,&Set[n].p,&Set[n].o);
            if (n % step == 0)
                printf("n = %i: (%i,%i,%i)\n",n,Set[n].s,Set[n].p,Set[n].o);
            if (Set[n].s > maxV) maxV = Set[n].s;
            if (Set[n].o > maxV) maxV = Set[n].o;
        }
        fclose(f);
        printf("Last one is %i,%i,%i\n",Set[n-1].s,Set[n-1].p,Set[n-1].o);
        printf("Maximum s,o = %i\n",maxV);
        printf("Writing...\n");
        for (j=0;j<S;j+=1000)
        { printf("j = %i...\n",j);
            for (i=1;i<=1000;i++)
            { char name[32];
                sprintf(name,"%04i.pairs",j+i);
                if (j+i < S) files[i] = fopen(name,"w");
            }
            for (n=0;n<n_triples;n++)
            { if (Set[n].p >= S)
                { printf("predicado %i fuera de rango\n",Set[n].p);
                    exit(1);
                }
                if ((Set[n].p > j) && (Set[n].p <= j+1000))
                { fwrite(&Set[n].s,sizeof(uint32_t),1,files[Set[n].p-j]);
                    fwrite(&Set[n].o,sizeof(uint32_t),1,files[Set[n].p-j]);
                }
            }
            for (i=1;i<=1000;i++)
                if (j+i < S) fclose(files[i]);
        }

}