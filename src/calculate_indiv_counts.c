//Xander Nuttle
//calculate_indiv_counts.c
//Calculates total number of reads for each individual from a .mipcounts file
//THIS PROGRAM IS NOT AT ALL GENERAL IN CURRENT FORM!!!
//Call: ./calculate_indiv_counts mipcounts_file

#include<stdio.h>
#include<string.h>

int main(int argc,char*argv[])
{
	FILE*countsfile;
	countsfile=fopen(*(argv+1),"r");
	char indiv[51];
	indiv[50]='\0';
	char newindiv[51];
	newindiv[50]='\0';
	char dummy[51];
	long coord;
	long numA,numB,numC,numD,numN;
	char ch;
	while((ch=getc(countsfile))!='\n')
	{
		continue;
	}
	long count=0;
	fpos_t pos;
	fgetpos(countsfile,&pos);
	fscanf(countsfile,"%s",indiv);
	fsetpos(countsfile,&pos);
	FILE*out;
	out=fopen("Troina2.indcts","w");
	while(fscanf(countsfile,"%s %s %ld %c %ld %ld %ld %ld %ld",newindiv,dummy,&coord,&ch,&numA,&numB,&numC,&numD,&numN)==9)
	{
		if(strcmp(indiv,newindiv)!=0)
		{
			fprintf(out,"%s\t%ld\n",indiv,count);
			count=0;
			strcpy(indiv,newindiv);
		}
		count+=(numA+numB+numC+numD+numN);
	}
	fprintf(out,"%s\t%ld\n",indiv,count);
	fclose(out);
	return 0;
}
