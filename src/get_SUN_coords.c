//Xander Nuttle
//get_SUN_coords.c
//Call: ./get_SUN_coords seqref.fasta seq2.fasta ... seqN.fasta alignment_base1_is_contig_base_N.bed gene_family.suns
//
//SUN characterization is important for identifying high-confidence single nucleotide markers that can be leveraged in experimental assays (with MIPs, for example)
//for paralog-specific genotyping. This program identifies SUNs from fasta files containing '-'s output from Jalview or another alignment editor and outputs their
//coordinates, both with respect to a master input fasta file and with respect to larger contig fasta files. The coordinate conversion is specified by a two column
//bed file, with the first column having larger contig names and the second the base in the corresponding larger contig that is the same base as the first base in the
//corresponding alignment output fasta file. All coordinates are in base 1.
//
//The first input file should be the alignment output fasta file containing the "master" sequence for the gene family from which any experimental assays (such as MIPs)
//will be designed -- since this is the ultimate purpose of SUN characterization, this is the only alignment output fasta coordinate system we need to use. However, we
//also need to use the larger contig coordinate systems to later identify SUNs intersecting unmasked SUNKs. Such SUNs can be assessed for presence/absence in several high
//coverage genomes sequenced with short-read technology. Thus, we will ultimately be able to identify high-quality SUNs present at high frequencies in normal human
//populations and link them with a given coordinate in the "master" sequence. These high-quality SUNs can then be given high priority during the experimental design process.
//The first column of the output file is the name of the contig having the corresponding SUN, the second column is the base 1 coordinate of the SUN with respect to its
//corresponding contig sequence, and the third column is the base 1 coordinate of the SUN with respect to the "master" sequence (from the alignment output, not from the
//"master" sequence's corresponding contig).


#include<stdio.h>
#include<string.h>
#include<stdlib.h>

int main(int argc,char*argv[])
{
	//open all files, read in bed file, store contig names and alignment start coordinates with respect to larger contig files
	FILE*files[argc-1];
	int i,j;
	for(i=0;i<(argc-2);i++)
	{
		files[i]=fopen(*(argv+1+i),"r");
	}
	files[i]=fopen(*(argv+1+i),"w");
	char contignames[i-1][51];
	char base[i];
	base[i-1]='\0';
	long alignstarts[i-1],aligncoords[i-1],delbases[i-1];
	j=i-1;
	for(i=0;i<j;i++)
	{
		contignames[i][50]='\0';
		delbases[i]=0;
	}
	i=0;
	while(fscanf(files[argc-3],"%s %ld",contignames[i],&(alignstarts[i]))==2)
		i++;

	//go through alignment output fasta files and output all SUNs: contig	larger_contig_coordinate	master_sequence_coordinate
	for(i=0;i<(argc-3);i++)
	{
		while((base[i]=getc(files[i]))!='\n')
			continue;
	}
	for(i=0;i<(argc-3);i++)
  {
    aligncoords[i]=1;
  }
	while((base[0]=getc(files[0]))!=EOF)
	{
		for(i=1;i<(argc-3);i++)
		{
			base[i]=getc(files[i]);
		}
		while(isspace(base[0]))
		{
			for(i=0;i<(argc-3);i++)
    	{
      	base[i]=getc(files[i]);
    	}
		}
		for(i=0;i<(argc-3);i++)
    {
      if(base[i]=='-')
			{
				delbases[i]++;
			}
    }
		for(i=0;i<(argc-3);i++)
    {
      if((strchr(base,base[i])==strrchr(base,base[i]))&&(strpbrk(base,"-YRWSKMDVHBN")==NULL))
      {
					fprintf(files[argc-2],"%s\t%ld\t%ld\n",contignames[i],aligncoords[i]-delbases[i]+alignstarts[i]-1,aligncoords[0]-delbases[0]);
			}
    }
		for(i=0;i<(argc-3);i++)
    {
        aligncoords[i]++;
    }
	}

	//cleanup and exit
	for(i=0;i<=(argc-2);i++)
	{
		fclose(files[i]);
	}
	return 0;
}
