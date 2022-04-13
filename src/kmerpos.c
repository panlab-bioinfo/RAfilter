#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include<math.h>
#include<stdint.h>
#include<time.h>
#include"uthash.h"

#define LINE_MAX 100000000
#define MASK 4398046511103


struct my_struct {
    long long int id;                    /* key */
    char name[30];
    UT_hash_handle hh;         /* makes this structure hashable */
};

struct my_struct *users = NULL;

void add_user(long long int *user_id, const char *name)
{
    struct my_struct *s;

    HASH_FIND_INT(users, user_id, s);  /* id already in the hash? */
    if (s == NULL) {
        s = (struct my_struct*)malloc(sizeof *s);
        s->id = *user_id;
        HASH_ADD_INT(users, id, s);  /* id is the key field */
    }
    strcpy(s->name, name);
}

struct my_struct *find_user(long long int *user_id)
{
    struct my_struct *s;

    HASH_FIND_INT(users, user_id, s);  /* s: output pointer */
    return s;
}

void delete_user(struct my_struct *user)
{
    HASH_DEL(users, user);  /* user: pointer to deletee */
    free(user);
}

void delete_all()
{
    struct my_struct *current_user;
    struct my_struct *tmp;

    HASH_ITER(hh, users, current_user, tmp) {
        HASH_DEL(users, current_user);  /* delete it (users advances to next) */
        free(current_user);             /* free it */
    }
}

int count_user() {
    return HASH_COUNT(users);
}

void print_users()
{
    struct my_struct *s;
	printf("size is %d\n", count_user(users));
    for (s = users; s != NULL; s = (struct my_struct*)(s->hh.next)) {
        printf("user id %llu: name %s\n", s->id, s->name);
    }
}

int basetonu(char c)
{
    switch (c)
        {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case'G':
            return 2;
		case 'T':
            return 3;
        }
}

unsigned long long int  kmertonu(char * a){
	int i;
	int v;
	unsigned long long int  k=0;
	char *b=(char*)malloc(sizeof(char)*strlen(a));
	for(i=0;i<strlen(a);i++)
	{
		v=basetonu(a[i]);
		k=v|(k<<2);
	}
    return k;
}

void read_kmer(char *path)
{
	FILE* fp = NULL;
	char buf[100]={0};
	unsigned long long int  num = 0;
	fp = fopen(path, "r");

	if (fp == NULL)
	{
		perror("my_fgets fopen");
		//return;
	}
	while (!feof(fp))
	{
		char* p = fgets(buf, sizeof(buf), fp);
		char *q=NULL; 
		if(p!=NULL){
		q = strtok(buf,"\t");
		num=kmertonu(q);
	    add_user(&num,q);
		}
		
	}
	if (fp != NULL)
	{
		fclose(fp);
		fp = NULL;
	}
}

char *getrc(char* kmer){
    char *kmer_rc=(char *)malloc(sizeof(char)*21);
	int i;
    for( i= 20;i >= 0;i--){
        char c = kmer[i];
        if(c == 'A' || c == 'a')
            kmer_rc[21 - i - 1] = 'T';
        else if(c == 'T' || c == 't')
            kmer_rc[21 - i - 1] = 'A';
        else if(c == 'C' || c == 'c')
            kmer_rc[21 - i - 1] = 'G';
        else if(c == 'G' || c == 'g')
            kmer_rc[21 - i - 1] = 'C';
        else{
           printf("ERROR!!!!!letters out of ATCG");
            exit(0);
        }
    }
    return kmer_rc;
}

void search_kmer(char *path)
{
	FILE* fp = NULL;//,*fw=NULL;
	char *buf=(char*)malloc(sizeof(char)*LINE_MAX);
	unsigned long long int  num = 0, b,c;
	char name[40]={0};
	fp = fopen(path, "r");
	//fw=fopen("1.txt","w+");
	if (fp == NULL)
	{
		perror("my_fgets fopen");
		return;
	}
	clock_t t1,t2,t3,t4,t5;
	while (!feof(fp))
	{
	 	
		char* p = fgets(buf, LINE_MAX, fp);
		int i=0,j=0,x=0;
		unsigned long long int f=0;
		int line_len =strlen(buf);
		if ('\n' == buf[line_len - 1]) {
			buf[line_len - 1] = '\0';
			line_len--;
		}
		char kmer[25]={0};
		struct my_struct *seq;
		if(p!=NULL){
		if(p[0]=='>')
		{
			strncpy(name,p,37);
			for(i=0;i<strlen(name);i++)
			{
				name[i]=name[i+1];
			}
			name[strlen(name)]='\0';
			
		}
		else
		{
			strncpy(kmer,p,21);
			long len=strlen(p);
			b=kmertonu(kmer);
			c=kmertonu(getrc(kmer));
			clock_t ts = clock();
			seq=find_user(&b);
		    if (seq!=NULL)
			{
				printf("%s\t%ld\t%d\t%s\t+\n",name,len,i,seq->name);
//				fprintf(fw,"%s\t%ld\t%d\t%s\t+\n",name,len,i,seq->name);
//				printf("%ld\n",ftell(fw));
			}
			seq=find_user(&c);
			if (seq!=NULL)
			{
				printf("%s\t%ld\t%d\t%s\t-\n",name,len,i,seq->name);
//				fprintf(fw,"%s\t%ld\t%d\t%s\t-\n",name,len,i,seq->name);
//				printf("%ld\n",ftell(fw));
			}
			printf("%ld",clock()-ts);
			i=21;
			while(i<len)
		    {
		        x=basetonu(p[i]);
		        // +forward
		        b=b<<2;
		        b=b&MASK;
		        b=b|x;
				seq=find_user(&b);
				clock_t st = clock();
		        if (seq!=NULL)
				{
					printf("%s\t%ld\t%d\t%s\t+\n",name,len,i-20,seq->name);
//					fprintf(fw,"%s\t%ld\t%d\t%s\t+\n",name,len,i-20,seq->name);
//					printf("%ld\n",ftell(fw));
				}
				printf("%lu",clock()-st);
				//-reverse
				c=c>>2;
				c=c&MASK;
		        f=~x;
		        f=f&3;
		        f=f<<40;
		        c=c|f;
		        seq=find_user(&c);
				if (seq!=NULL)
				{
					printf("%s\t%ld\t%d\t%s\t-\n",name,len,i-20,seq->name);
//					fprintf(fw,"%s\t%ld\t%d\t%s\t-\n",name,len,i-20,seq->name);
//					printf("%ld\n",ftell(fw));
				}

		        i++;
		    }   
    
		}
		}
	}
	//fclose(fw);
	if (fp != NULL)
	{
		fclose(fp);
		fp = NULL;
	}
}

int main(int argc,char* argv[])
{
	char *kmer="/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.kmers.1.0_100M.dump";
    	char *seq="/data/panweihua/human_chm13/hifi/fastqs/chm13.asm.p_ctg.compressed.fa";
//    char *kmer=argv[1];
//    char *seq=argv[2];
//    char *kmer="unique_kmer.txt";//argv[1];
//    char *seq="reads.txt";
	read_kmer(kmer);
	search_kmer(seq);
	delete_all();
	return 0;
}
