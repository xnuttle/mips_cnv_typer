/*
 * Xander Nuttle
 *
 * call_mip_srgap2_cn.c
 *
 * Generates paralog-specific SRGAP2 copy number calls from MIP paralog-specific
 * read count data, using a maximum likelihood approach and dynamic programming.
 *
 * Call: ./call_mip_srgap2_cn SRGAP2_miptargets_file SRGAP2_mip_read_counts_file output_file_base_name
 * Example call: ./call_mip_srgap2_cn SRGAP2_final.miptargets pos_ctrl_expt_SRGAP2_final.mipcounts pos_ctrl_expt
 *
 * NOTE: The mipcounts file must contain data for only MIPs in the miptargets
 * file!
 *
 * This program takes in (1) a file containing information on SRGAP2 MIP targets
 * and (2) a file containing paralog-specific SRGAP2 read counts, outputting (1)
 * a file containing SRGAP2 copy number genotypes for each individual included
 * in the paralog-specific SRGAP2 read counts file. Input file (2) should
 * contain data corresponding only to MIPs in input file (1); that is, the data
 * in file 2 should be generated by running mip_analysis.c with file 1 as an
 * input.  Additionally, all MIPs in file (1) must be sorted according to master
 * sequence target coordinate, smallest to largest.
 *
 * This program generates genotype calls using a maximum likelihood approach
 * with dynamic programming. For each individual, likelihoods of observing the
 * count data for each MIP are calculated under 400 different possible
 * underlying SRGAP2 paralog-specific copy number states, where SRGAP2A and
 * SRGAP2C can have copy number 0-3 and SRGAP2B and SRGAP2D can have copy number
 * 0-4 (400 combinations = 4*5*4*5). Then, for each individual, these
 * likelihoods are used to construct a weighted directed acyclic graph (WDAG),
 * with prior probabilities based on observed genotype data from previous
 * experiments incorporated into the likelihoods for the first MIP count
 * data. As the graph is constructed, highest scoring paths allowing 1 and 2
 * transitions between copy number states are tracked. Transitions are
 * restricted to duplications or deletions affecting 1 paralog only or gene
 * conversion events, where the copy numbers of two paralogs change, but the
 * total number of SRGAP2 copies remains constant. Placing these restrictions on
 * transitions reflects the fact that true biological events should fall into
 * one of these categories (single-paralog-affecting duplicaion/deletion or gene
 * conversion). The program ultimately identifies the highest-scoring (most
 * likely) paths through the likelihood graph allowing 0, 1, and 2 transitions
 * and their corresponding scores (likelihoods). A heuristic is used to assess
 * increases in likelihood of the 1-transition and 2-transition paths compared
 * to the 0 transition path and determine whether they signal an internal event
 * (duplication, deletion, or conversion) and warrant calling the genotype for
 * an individual as complex. In most cases, the scores of the 1-transition and
 * 2-transition paths will not be substantially higher than that of the
 * 0-transition path, and the genotype for an individual will be called as
 * simple (a single copy number state across the entirety of the gene).
 */

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "iniparser.h"

/*
 * Define constants.
 */
#define PARALOG_COUNTS_SECTION "paralog_copy_number_states"
#define MAX_STRING_LENGTH 80
#define BUFFER_SIZE 1000

// Parameter setting the minimum likelihood value for a single data point
#define MIN_LIKELIHOOD -30

static const long SRGAP2D_DEL_START=105947; //master sequence base 1 coordinate of SRGAP2D deletion region start
static const long SRGAP2D_DEL_END=213356; //master sequence base 1 coordinate of SRGAP2D deletion region end
static const int MIN_LIKELIHOOD_DIFF=40; //see below for detailed description of this heuristic
static const int MIN_MIPS_IN_CN_STATE=5; //see below for detailed description of this heuristic
static const int num_args=4; //number of required command line arguments

/*
 * Given a pointer to an iniparser instance, load the paralog copy number states
 * for all paralogs into a vector.
 */
gsl_vector_int* get_paralog_copy_numbers(dictionary* ini) {
    gsl_vector_int* paralog_copy_numbers = NULL;
    int total_paralogs, i;
    char** paralog_keys;

    /*
     * Get total number of paralogs defined in the configuration file and keys
     * for each paralog entry in the file.
     */
    total_paralogs = iniparser_getsecnkeys(ini, PARALOG_COUNTS_SECTION);
    paralog_keys = iniparser_getseckeys(ini, PARALOG_COUNTS_SECTION);
    printf("Found %i paralogs\n", total_paralogs);

    if (total_paralogs == 0) {
        return NULL;
    }

    /*
     * Allocate a vector for the paralog copy numbers and populate the vector
     * from the configuration file using the configuration keys.
     */
    paralog_copy_numbers = gsl_vector_int_alloc(total_paralogs);
    for (i = 0; i < total_paralogs; i++) {
        gsl_vector_int_set(paralog_copy_numbers, i, iniparser_getint(ini, paralog_keys[i], -1));
    }

    return paralog_copy_numbers;
}

/*
 * Given a pointer to an iniparser instance, load priors for the given paralog
 * number.
 */
gsl_vector* get_paralog_priors(dictionary* ini, int paralog, int paralog_copy_states) {
    char section[MAX_STRING_LENGTH];
    char** paralog_keys;
    double total, prior;
    int i;
    gsl_vector* priors;

    sprintf(section, "priors_paralog%i", paralog);

    paralog_keys = iniparser_getseckeys(ini, section);
    if (paralog_keys == NULL) {
        return NULL;
    }

    total = iniparser_getdouble(ini, paralog_keys[0], -1);
    priors = gsl_vector_alloc(paralog_copy_states);

    for (i = 0; i < paralog_copy_states; i++) {
        prior = iniparser_getdouble(ini, paralog_keys[i + 1], -1);

        // Don't let prior values equal 0 for statistical reasons.
        if (prior == 0) {
            prior = exp(MIN_LIKELIHOOD);
        }
        else {
            prior = prior / total;
        }

        gsl_vector_set(priors, i, prior);
    }

    return priors;
}

/*
 * Given n values and an array x of integers for each n value, populate a matrix
 * with n columns and rows equal to the product of the length of the n-1 matrix
 * and x[n].
 */
gsl_matrix_int* populate_matrix(int n, gsl_vector_int* x) {
    int i, j, k;
    gsl_matrix_int* matrix;

    if (n == 1) {
        matrix = gsl_matrix_int_alloc(gsl_vector_int_get(x, n - 1), n);

        for (i=0; i < gsl_vector_int_get(x, n - 1); i++) {
            gsl_matrix_int_set(matrix, i, n-1, i);
        }
    }
    else {
        gsl_matrix_int* previous_matrix = populate_matrix(n-1, x);
        int previous_rows = (int)previous_matrix->size1;
        matrix = gsl_matrix_int_alloc(previous_rows*gsl_vector_int_get(x, n - 1), n);

        /*
         * Fill in the current matrix with i copies of the previous matrix's
         * values which fill up to n - 2. Populate the final column with values
         * for this paralog.
         */
        for (i = 0; i < gsl_vector_int_get(x, n - 1); i++) {
            for (j = 0; j < previous_rows; j++) {
                for(k = 0; k < n - 1; k++) {
                    gsl_matrix_int_set(
                        matrix,
                        previous_rows * i + j,
                        k,
                        gsl_matrix_int_get(previous_matrix, j, k)
                    );
                }

                // Populate the final column with values for this paralog.
                gsl_matrix_int_set(matrix, previous_rows * i + j, n-1, i);
            }
        }

        gsl_matrix_int_free(previous_matrix);
    }

    return matrix;
}

int main(int argc,char*argv[])
{
    long i, j;

    // Check to make sure there are there are enough command line arguments provided.
    if(argc<(num_args+1)) //argc includes the program call in its count
    {
        printf("Usage: %s miptargets_file mip_read_counts_file output_file_base_name configuration_file\n\n", argv[0]);
        printf("Example call: %s SRGAP2_final.miptargets pos_ctrl_expt_SRGAP2_final.mipcounts pos_ctrl_expt SRGAP2.conf.ini\n", argv[0]);
        return 1;
    }

    /*
     * Load the configuration file.
     */
    dictionary* ini;
    ini = iniparser_load(argv[4]);
    if (ini == NULL) {
        fprintf(stderr, "Cannot open configuration file: %s\n", argv[4]);
        return -1;
    }

    // Load paralog copy numbers.
    gsl_vector_int* paralog_copy_numbers = get_paralog_copy_numbers(ini);
    if (paralog_copy_numbers == NULL) {
        fprintf(stderr, "Couldn't load paralog copy numbers.\n");
        iniparser_freedict(ini);
        return -1;
    }

    int number_of_paralogs = (int)paralog_copy_numbers->size;
    int number_of_copy_states = 1;
    for (i = 0; i < number_of_paralogs; i++) {
        number_of_copy_states = number_of_copy_states * gsl_vector_int_get(paralog_copy_numbers, i);
    }

    /*
     * Load paralog priors. Priors give the probabilities of observing a copy
     * number state of 0,1,2,3...n. for each paralog. These probabilities are
     * based on empirical data from screening and genotyping efforts.
     *
     * Extremely unlikely states (never observed so far but hypothetically
     * possible) are arbitrarily assigned a value with a resulting natural
     * logarithm of -30.
     */
    gsl_vector* priors_by_paralog[number_of_paralogs];
    for (i = 0; i < number_of_paralogs; i++) {
        priors_by_paralog[i] = get_paralog_priors(ini, i, gsl_vector_int_get(paralog_copy_numbers, i));
        if (priors_by_paralog[i] == NULL) {
            fprintf(stderr, "Couldn't load priors for paralog %ld.\n", i);
            return -1;
        }
    }

    // Clean up configuration file.
    iniparser_freedict(ini);

    // Get information about number of mip targets designed for copy number genotyping.
    FILE*miptargetsfile;
    fpos_t pos;
    char miptype;
    char dummy[501];
    long num_mip_targets=0,num_exon_mips=0;
    miptargetsfile=fopen(*(argv+1),"r");
    fgetpos(miptargetsfile,&pos);
    while(fscanf(miptargetsfile,"%s %s %s %c %s %s %s %s %s",dummy,dummy,dummy,&miptype,dummy,dummy,dummy,dummy,dummy)==9)
    {
        if(miptype!='E')
        {
            num_mip_targets++;
        }
        else
        {
            num_exon_mips++;
        }
    }
    fsetpos(miptargetsfile,&pos);

    /*
     * Get and store information about the gene family's master sequence target
     * coordinate and specificity for each MIP.
     */
    char spec_string[number_of_paralogs];
    // A string to hold either 0 or 1 (a MIP either has or lacks specificity for a given paralog).
    char spec[2];
    long start,end;
    long target_coords[num_mip_targets];
    long specificities[num_mip_targets][number_of_paralogs];
    i=0;
    spec[1]='\0';
    while(fscanf(miptargetsfile,"%s %ld %c %ld %s %c %s %s %s %s %s",dummy,&start,&miptype,&end,dummy,&miptype,dummy,dummy,spec_string,dummy,dummy)==11)
    {
        if(miptype!='E')
        {
            target_coords[i]=(start+end)/2;
            for(j=0;j<number_of_paralogs;j++)
            {
                spec[0]=spec_string[j];
                specificities[i][j]=strtol(spec,NULL,10);
            }
            i++;
        }
    }
    fclose(miptargetsfile);

    /*
     * Generate vector of possible copy number states and vector of ln(prior
     * probabilities) corresponding to each copy number state assume
     * independence of paralog-specific copy number genotypes
     */
    int copy_states[number_of_copy_states][number_of_paralogs];
    // Vector of prior probabilities for each copy number state.
    double priors[number_of_copy_states];

    /*
     * Populate copy number states into a matrix to allow dynamic configuration
     * of paralogs. Load results into the copy states array to maintain a
     * consistent interface in the following code.
     */
    gsl_matrix_int* matrix = populate_matrix(number_of_paralogs, paralog_copy_numbers);
    for (i = 0; i < (int)matrix->size1; i++) {
        for (j = 0; j < (int)matrix->size2; j++) {
            copy_states[i][j] = gsl_matrix_int_get(matrix, i, j);
        }
    }

    /*
     * Free memory allocated for the matrix now that its contents have been
     * copied into the copy states array.
     */
    if (matrix != NULL) {
        gsl_matrix_int_free(matrix);
    }

    /*
     * Initialize prior probabilities.
     */
    double prior;
    for (i = 0; i < number_of_copy_states; i++) {
        prior = 0;
        for (j = 0; j < number_of_paralogs; j++) {
            prior = prior + log(gsl_vector_get(priors_by_paralog[j], copy_states[i][j]));
        }

        priors[i] = prior;
    }

    // Free memory allocated for paralog prior vectors.
    for (i = 0; i < number_of_paralogs; i++) {
        if (priors_by_paralog[i] != NULL) {
            gsl_vector_free(priors_by_paralog[i]);
        }
    }

    /*
     * Set up node structure for dynamic programming to find maximum likelihood
     * path through graph transitions to allow = 2*(number of internal events to
     * allow), because each internal event can maximally contribute two edges in
     * copy number state corresponding to 2 transitions right now, the program
     * allows for 1 internal event (testing to see if this will still pick up >1
     * internal event, in which case allowing for more internal events will not
     * be necessary)
     */
    struct node
    {
        long state;
        double max0;
        double max1;
        double max2;
        struct node*prev_node_1trans;
        struct node*prev_node_2trans;
        long index;
    };

    // Allocate storage for likelihood graph.
    struct node*likelihood_graph;
    likelihood_graph=(struct node*)malloc(num_mip_targets*number_of_copy_states*sizeof(struct node));

    // Initialize output files.
    FILE*out,*out2,*out3;
    char output_base_name[51];
    char output_file_extension[9]=".cncalls";
    char ext2[12]=".compevents";
    char ext3[13]=".simplecalls";
    strncpy(output_base_name,*(argv+3),42);
    strcat(output_base_name,output_file_extension);
    out=fopen(output_base_name,"w");
    for(i=0;i<51;i++)
    {
        output_base_name[i]='\0';
    }
    strncpy(output_base_name,*(argv+3),39);
    strcat(output_base_name,ext2);
    out2=fopen(output_base_name,"w");
    for(i=0;i<51;i++)
    {
        output_base_name[i]='\0';
    }
    strncpy(output_base_name,*(argv+3),38);
    strcat(output_base_name,ext3);
    out3=fopen(output_base_name,"w");

    fprintf(out3,"Individual");
    for (i = 0; i < number_of_paralogs; i++) {
        fprintf(out3, "\tParalog%ld_CN", i);
    }
    fprintf(out3, "\tPossible_Complex_CN_Genotype\n");

    // For each individual, read in counts, calculate and store individual likelihoods of data for each MIP under each copy number state.
    FILE*countsfile;
    char individual[31];
    long counts[number_of_paralogs+1];
    long indiv_counts[num_mip_targets][number_of_paralogs+1];
    double probs[number_of_paralogs+1];
    double L;
    long k, p;
    double mip_likelihoods[num_mip_targets][number_of_copy_states];

    char* delimiter = "\t";
    char* token = NULL;
    char* strtol_end;
    char buffer[BUFFER_SIZE];
    char* token_buffer;

    individual[30]='\0';
    countsfile=fopen(*(argv+2),"r");
    while(getc(countsfile)!='\n') {
        continue;
    }
    fgetpos(countsfile,&pos);

    /*
     * Parse rows for each sample in the counts file. Each sample is processed
     * in one iteration of this loop.
     */
    while(fgets(buffer, sizeof buffer, countsfile))
    {
        fsetpos(countsfile,&pos);
        j=0;
        for(i=0;i<(num_mip_targets+num_exon_mips);i++)
        {
            // Parse this line from the counts file.
            fgets(buffer, sizeof buffer, countsfile);

            // First column is the sample.
            token = strtok_r(buffer, delimiter, &token_buffer);
            strcpy(individual, token);

            // Second column is the target sequence.
            token = strtok_r(NULL, delimiter, &token_buffer);
            strcpy(dummy, token);

            // Third column is the MIP coordinate which can be skipped here..
            token = strtok_r(NULL, delimiter, &token_buffer);

            // Fourth column is the MIP type.
            token = strtok_r(NULL, delimiter, &token_buffer);
            miptype = token[0];

            // All remaining columns should be counts for each paralog plus an "uncertainty" bin.
            for (k = 0; k < number_of_paralogs + 1; k++) {
                token = strtok_r(NULL, delimiter, &token_buffer);
                counts[k] = strtol(token, &strtol_end, 10);
            }

            if(miptype!='E')
            {
                for(k=0;k<(number_of_paralogs+1);k++)
                {
                    indiv_counts[j][k]=counts[k];
                }
                j++;
            }
        }

        // Calculate likelihoods of observed counts at each MIP under each possible copy number state.
        double sum_of_copy_states;
        for(i=0;i<num_mip_targets;i++)
        {
            for(j=0;j<number_of_copy_states;j++)
            {
                sum_of_copy_states = 0;
                for (k = 0; k < number_of_paralogs; k++) {
                    sum_of_copy_states = sum_of_copy_states + (double)copy_states[j][k];
                }

                probs[number_of_paralogs]=0.0;
                // MIP not in SRGAP2D deletion region
                // TODO: remove this SRGAP2-specific code.
                if((target_coords[i]<SRGAP2D_DEL_START)||(target_coords[i]>SRGAP2D_DEL_END))
                {
                    for(k=0;k<number_of_paralogs;k++)
                    {
                        if(specificities[i][k]==1)
                        {
                            probs[k]=(double)(copy_states[j][k]) / sum_of_copy_states;
                        }
                        else
                        {
                            probs[k]=0.0;
                            probs[number_of_paralogs]+=(double)(copy_states[j][k]) / sum_of_copy_states;
                        }
                    }
                }
                else
                {
                    sum_of_copy_states = 0;
                    for (k = 0; k < number_of_paralogs - 1; k++) {
                        sum_of_copy_states = sum_of_copy_states + (double)copy_states[j][k];
                    }

                    probs[number_of_paralogs-1]=0.0;
                    for(k=0;k<(number_of_paralogs-1);k++)
                    {
                        if(specificities[i][k]==1)
                        {
                            probs[k]=(double)(copy_states[j][k]) / sum_of_copy_states;
                        }
                        else
                        {
                            probs[k]=0.0;
                            probs[number_of_paralogs]+=(double)(copy_states[j][k]) / sum_of_copy_states;
                        }
                    }
                }

                // Initialize countvec.
                unsigned int countvec[number_of_paralogs + 1];
                for (p = 0; p < number_of_paralogs + 1; p++) {
                    countvec[p] = indiv_counts[i][p];
                }

                L=gsl_ran_multinomial_lnpdf(number_of_paralogs + 1,probs,countvec);
                if(L<MIN_LIKELIHOOD)
                {
                    // Minimum log likelihood value for 1 MIP probe arbitrarily assigned to -30.
                    L=MIN_LIKELIHOOD;
                }
                // i = mip_target, j = copy_state
                mip_likelihoods[i][j]=L;
            }
        }

        /*
         * Use dynamic programming to determine whether the data support more
         * than 1 copy number state across the gene family.
         */
        double max_max0=-DBL_MAX,max_max1=-DBL_MAX,max_max2=-DBL_MAX;
        long maxindex_max0,maxindex_max1,maxindex_max2;

        // Setup first column and calculate maximum likelihood value and index of corresponding node.
        for(j=0;j<number_of_copy_states;j++)
        {
            likelihood_graph[j].index=j;
            likelihood_graph[j].state=j;
            likelihood_graph[j].max0=mip_likelihoods[0][j]+priors[j];
            likelihood_graph[j].max1=-DBL_MAX;
            likelihood_graph[j].max2=-DBL_MAX;
            likelihood_graph[j].prev_node_1trans=NULL;
            likelihood_graph[j].prev_node_2trans=NULL;
        }

        // Setup all other verticies.
        for(j=number_of_copy_states;j<(number_of_copy_states*num_mip_targets);j++)
        {
            if((j%number_of_copy_states)==0)
            {
                max_max0=-DBL_MAX;
                max_max1=-DBL_MAX;
                for((k=j-number_of_copy_states);k<j;k++)
                {
                    if(likelihood_graph[k].max0>max_max0)
                    {
                        max_max0=likelihood_graph[k].max0;
                        maxindex_max0=k;
                    }
                    if(likelihood_graph[k].max1>max_max1)
                    {
                        max_max1=likelihood_graph[k].max1;
                        maxindex_max1=k;
                    }
                }
            }
            likelihood_graph[j].index=j;
            likelihood_graph[j].state=j%number_of_copy_states;
            likelihood_graph[j].max0=mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states]+likelihood_graph[j-number_of_copy_states].max0;
            if(max_max0>likelihood_graph[j-number_of_copy_states].max1)
            {
                likelihood_graph[j].max1=max_max0+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                likelihood_graph[j].prev_node_1trans=&(likelihood_graph[maxindex_max0]);
            }
            else
            {
                likelihood_graph[j].max1=likelihood_graph[j-number_of_copy_states].max1+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                likelihood_graph[j].prev_node_1trans=&(likelihood_graph[j-number_of_copy_states]);
            }
            if(max_max1>likelihood_graph[j-number_of_copy_states].max2)
            {
                likelihood_graph[j].max2=max_max1+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                likelihood_graph[j].prev_node_2trans=&(likelihood_graph[maxindex_max1]);
            }
            else
            {
                likelihood_graph[j].max2=likelihood_graph[j-number_of_copy_states].max2+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                likelihood_graph[j].prev_node_2trans=&(likelihood_graph[j-number_of_copy_states]);
            }
        }

        /*
         * Determine maximum values of final likelihoods allowing 0, 1, and 2
         * transitions and use to determine whether there may be multiple
         * underlying copy number states across the gene family.
         */
        max_max0=-DBL_MAX;
        max_max1=-DBL_MAX;
        max_max2=-DBL_MAX;
        for(j=(number_of_copy_states*num_mip_targets-number_of_copy_states);j<(number_of_copy_states*num_mip_targets);j++)
        {
            if(likelihood_graph[j].max0>max_max0)
            {
                max_max0=likelihood_graph[j].max0;
                // Maximally likely copy number state assuming no internal events.
                maxindex_max0=j%number_of_copy_states;
            }
            if(likelihood_graph[j].max1>max_max1)
            {
                max_max1=likelihood_graph[j].max1;
                maxindex_max1=j;
            }
            if(likelihood_graph[j].max2>max_max2)
            {
                max_max2=likelihood_graph[j].max2;
                maxindex_max2=j;
            }
        }
        fprintf(out,"Individual: %s\nL0=%lf L1-L0=%lf L2-L0=%lf\n",individual,max_max0,max_max1-max_max0,max_max2-max_max0);
        fprintf(out,"Genotype (assuming no internal events): ");
        for (i = 0; i < number_of_paralogs; i++) {
            fprintf(out, "%d", copy_states[maxindex_max0][i]);
        }
        fprintf(out, "\n");

        int evidence_for_1T=0,evidence_for_2T=0;
        long num_mips_state1=0,num_mips_state2=0,num_mips_state3=0;
        struct node*current;
        long oldstate;
        int trans_made=0;

        /*
         * Heuristics for predicting 1 copy number state transition
         *  Likelihood_1_transition - Likelihood_0_transitions > 40 (MIN_LIKELIHOOD_DIFF)
         *  At least 5 MIPs in each copy number state (MIN_MIPS_IN_CN_STATE)
         *  Likelihood_2_transitions - Likelihood_1_transition <= 40
         *
         * Heuristics for predicting 2 copy number state transitions
         *  Likelihood_2_transitions - Likelihood_0_transitions > 40
         *  At least 5 MIPs in each copy number state
         *  Likelihood_2_transitions - Likelihood_1_transition > 40
         */
        if((max_max1-max_max0)>MIN_LIKELIHOOD_DIFF)
        {
            // Trace back to determine number of MIPs in each state.
            // 1 transition
            current=&(likelihood_graph[maxindex_max1]);
            while(1)
            {
                num_mips_state2++;
                oldstate=current->state;
                current=current->prev_node_1trans;
                if((current->state)!=oldstate)
                {
                    break;
                }
            }
            num_mips_state1=num_mip_targets-num_mips_state2;
            if((num_mips_state1>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state2>=MIN_MIPS_IN_CN_STATE))
            {
                evidence_for_1T=1;
            }
        }
        if((max_max2-max_max0)>MIN_LIKELIHOOD_DIFF)
        {
            // 2 transitions
            num_mips_state1=0;
            num_mips_state2=0;
            current=&(likelihood_graph[maxindex_max2]);
            while(trans_made<2)
            {
                if((max_max2-max_max1)==0)
                {
                    break;
                }

                oldstate=current->state;
                if(trans_made==0)
                {
                    num_mips_state3++;
                    current=current->prev_node_2trans;
                }
                else
                {
                    num_mips_state2++;
                    current=current->prev_node_1trans;
                }
                if((current->state)!=oldstate)
                {
                    trans_made++;
                }
            }
            num_mips_state1=num_mip_targets-num_mips_state2-num_mips_state3;
            if((num_mips_state1>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state2>=MIN_MIPS_IN_CN_STATE)&&(num_mips_state3>=MIN_MIPS_IN_CN_STATE))
            {
                if(evidence_for_1T)
                {
                    if(max_max2-max_max1>MIN_LIKELIHOOD_DIFF)
                    {
                        evidence_for_2T=1;
                        evidence_for_1T=0;
                    }
                }
                else
                {
                    evidence_for_2T=1;
                    evidence_for_1T=0;
                }
            }
        }

        /*
         * If there is evidence for possible multiple copy number states across
         * the gene family, use dynamic programming to find maximum likelihood
         * paths through likelihood graph, allowing up to 2 transitions (1
         * internal event).
         */
        int num_paralog_cn_changes,total_copies_prev,total_copies;
        if(evidence_for_1T||evidence_for_2T)
        {
            // Alert user.
            printf("\n%s may have a complex genotype!!! Graphically view the data for %s to ensure genotyping accuracy.\n",individual,individual);
            // Print to 2nd output file so you can easily visualize only these individuals.
            fprintf(out2,"%s\t",individual);
            fprintf(out2,"L0=%lf L1-L0=%lf L2-L0=%lf\t",max_max0,max_max1-max_max0,max_max2-max_max0);

            fprintf(out2, "Genotype (assuming no internal events): ");
            fprintf(out3, "%s", individual);
            for (i = 0; i < number_of_paralogs; i++) {
                fprintf(out2, "%d", copy_states[maxindex_max0][i]);
                fprintf(out3, "\t%d", copy_states[maxindex_max0][i]);
            }
            fprintf(out2, "\t");
            fprintf(out3, "\tYES\n");

            // Setup first column and calculate maximum likelihood value and index of corresponding node.
            for(j=0;j<number_of_copy_states;j++)
            {
                likelihood_graph[j].state=j;
                likelihood_graph[j].index=j;
                likelihood_graph[j].max0=mip_likelihoods[0][j]+priors[j];
                likelihood_graph[j].max1=-DBL_MAX;
                likelihood_graph[j].max2=-DBL_MAX;
                likelihood_graph[j].prev_node_1trans=NULL;
                likelihood_graph[j].prev_node_2trans=NULL;
            }
            // Setup all other verticies.
            for(j=number_of_copy_states;j<(number_of_copy_states*num_mip_targets);j++)
            {
                max_max0=-DBL_MAX;
                max_max1=-DBL_MAX;
                maxindex_max0=0;
                maxindex_max1=0;
                for((k=j-number_of_copy_states);k<j;k++)
                {
                    num_paralog_cn_changes=0;
                    total_copies_prev=0;
                    total_copies=0;
                    for(i=0;i<4;i++)
                    {
                        if(copy_states[k%number_of_copy_states][i]!=copy_states[j%number_of_copy_states][i])
                        {
                            num_paralog_cn_changes++;
                            total_copies_prev+=copy_states[k%number_of_copy_states][i];
                            total_copies+=copy_states[j%number_of_copy_states][i];
                        }
                    }

                    // Signature of gene conversion, still a single event.
                    if((num_paralog_cn_changes==2)&&(total_copies_prev==total_copies))
                    {
                        num_paralog_cn_changes=1;
                    }

                    // Only allow a single event.
                    if((likelihood_graph[k].max0>max_max0)&&(num_paralog_cn_changes==1))
                    {
                        max_max0=likelihood_graph[k].max0;
                        maxindex_max0=k;
                    }

                    // Only allow a single event.
                    if((likelihood_graph[k].max1>max_max1)&&(num_paralog_cn_changes==1))
                    {
                        max_max1=likelihood_graph[k].max1;
                        maxindex_max1=k;
                    }
                }
                likelihood_graph[j].state=j%number_of_copy_states;
                likelihood_graph[j].index=j;
                likelihood_graph[j].max0=mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states]+likelihood_graph[j-number_of_copy_states].max0;
                if(max_max0>likelihood_graph[j-number_of_copy_states].max1)
                {
                    likelihood_graph[j].max1=max_max0+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                    likelihood_graph[j].prev_node_1trans=&(likelihood_graph[maxindex_max0]);
                }
                else
                {
                    likelihood_graph[j].max1=likelihood_graph[j-number_of_copy_states].max1+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                    likelihood_graph[j].prev_node_1trans=&(likelihood_graph[j-number_of_copy_states]);
                }
                if(max_max1>likelihood_graph[j-number_of_copy_states].max2)
                {
                    likelihood_graph[j].max2=max_max1+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                    likelihood_graph[j].prev_node_2trans=&(likelihood_graph[maxindex_max1]);
                }
                else
                {
                    likelihood_graph[j].max2=likelihood_graph[j-number_of_copy_states].max2+mip_likelihoods[j/number_of_copy_states][j%number_of_copy_states];
                    likelihood_graph[j].prev_node_2trans=&(likelihood_graph[j-number_of_copy_states]);
                }
            }

            // Determine maximum values of final likelihoods allowing 0, 1, and 2 transitions.
            max_max0=-DBL_MAX;
            max_max1=-DBL_MAX;
            max_max2=-DBL_MAX;
            maxindex_max2=0;
            for(j=(number_of_copy_states*num_mip_targets-number_of_copy_states);j<(number_of_copy_states*num_mip_targets);j++)
            {
                if(likelihood_graph[j].max0>max_max0)
                {
                    max_max0=likelihood_graph[j].max0;
                    // Maximally likely copy number state assuming no internal events.
                    maxindex_max0=j%number_of_copy_states;
                }
                if(likelihood_graph[j].max1>max_max1)
                {
                    max_max1=likelihood_graph[j].max1;
                    maxindex_max1=j;
                }
                if(likelihood_graph[j].max2>max_max2)
                {
                    max_max2=likelihood_graph[j].max2;
                    maxindex_max2=j;
                }
            }

            // Trace back to determine copy number states at each MIP for 1 and 2 transition scenarios.
            struct node*current;
            // 1 transition
            long first_mip_state2_1trans,state1_1trans,oldstate,oldindex;
            current=&(likelihood_graph[maxindex_max1]);
            while((current->index)>=number_of_copy_states)
            {
                oldstate=current->state;
                oldindex=current->index;
                current=current->prev_node_1trans;
                if((current->state)!=oldstate)
                {
                    first_mip_state2_1trans=target_coords[oldindex/number_of_copy_states];
                    state1_1trans=current->state;
                    break;
                }
            }

            // 2 transitions
            long first_mip_state2_2trans,last_mip_state2_2trans,state2_2trans,state1_2trans;
            current=&(likelihood_graph[maxindex_max2]);
            int trans_made=0;
            while(trans_made<2)
            {
                oldstate=current->state;
                oldindex=current->index;
                if(trans_made==0)
                {
                    current=current->prev_node_2trans;
                }
                else
                {
                    current=current->prev_node_1trans;
                }
                if((current->state)!=oldstate)
                {
                    if(trans_made==0)
                    {
                        state2_2trans=current->state;
                        last_mip_state2_2trans=target_coords[(current->index)/number_of_copy_states];
                        trans_made++;
                    }
                    else
                    {
                        state1_2trans=current->state;
                        first_mip_state2_2trans=target_coords[oldindex/number_of_copy_states];
                        trans_made++;
                    }
                }
            }

            if(evidence_for_1T)
            {
                fprintf(out,"Genotype (assuming 1 transition): ");
                fprintf(out2,"Genotype (assuming 1 transition): ");

                for (i = 0; i < number_of_paralogs; i++) {
                    fprintf(out, "%d", copy_states[state1_1trans][i]);
                    fprintf(out2, "%d", copy_states[state1_1trans][i]);
                }

                fprintf(out, " ");
                fprintf(out2, " ");

                for (i = 0; i < number_of_paralogs; i++) {
                    fprintf(out, "%d", copy_states[maxindex_max1 % number_of_copy_states][i]);
                    fprintf(out2, "%d", copy_states[maxindex_max1 % number_of_copy_states][i]);
                }

                fprintf(out, ", First MIP in State 2: %ld\n\n", first_mip_state2_1trans);
                fprintf(out2, ", First MIP in State 2: %ld\n\n", first_mip_state2_1trans);
            }
            else
            {
                fprintf(out, "Genotype (assuming 2 transitions): ");
                fprintf(out2, "Genotype (assuming 2 transitions): ");

                for (i = 0; i < number_of_paralogs; i++) {
                    fprintf(out, "%d", copy_states[state1_2trans][i]);
                    fprintf(out2, "%d", copy_states[state1_2trans][i]);
                }

                fprintf(out, " ");
                fprintf(out2, " ");

                for (i = 0; i < number_of_paralogs; i++) {
                    fprintf(out, "%d", copy_states[state2_2trans][i]);
                    fprintf(out2, "%d", copy_states[state2_2trans][i]);
                }

                fprintf(out, " ");
                fprintf(out2, " ");

                for (i = 0; i < number_of_paralogs; i++) {
                    fprintf(out, "%d", copy_states[maxindex_max2 % number_of_copy_states][i]);
                    fprintf(out2, "%d", copy_states[maxindex_max2 % number_of_copy_states][i]);
                }

                fprintf(out, ", First MIP in State 2: %ld, Last MIP in State 2: %ld\n\n", first_mip_state2_2trans, last_mip_state2_2trans);
                fprintf(out2, ", First MIP in State 2: %ld, Last MIP in State 2: %ld\n\n", first_mip_state2_2trans, last_mip_state2_2trans);
            }
        }
        else
        {
            fprintf(out,"\n");
            fprintf(out3, "%s", individual);
            for (i = 0; i < number_of_paralogs; i++) {
                fprintf(out3, "\t%d", copy_states[maxindex_max0][i]);
            }
            fprintf(out3, "\tNO\n");
        }

        // Prepare for next individual.
        fgetpos(countsfile,&pos);
    }

    // Clean up and exit.
    free(likelihood_graph);
    fclose(out);
    fclose(out2);
    fclose(out3);
    return 0;
}
