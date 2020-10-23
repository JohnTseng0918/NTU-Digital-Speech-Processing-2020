#include "../inc/hmm.h"
#include <stdio.h>

double alpha[MAX_SEQ][MAX_STATE];
double beta[MAX_SEQ][MAX_STATE];
double Gamma[MAX_SEQ][MAX_STATE];
double epsilon[MAX_STATE][MAX_STATE];

double initial_sum[MAX_STATE];
double epsilon_sum[MAX_STATE][MAX_STATE];
double epsilon_sum_div[MAX_STATE][MAX_STATE];
double Gamma_sum[MAX_SEQ][MAX_STATE];
double Gamma_sum_mol[MAX_SEQ][MAX_STATE];

int seq_to_num(char c)
{
    return c - 'A';
}

void init_table()
{
    int i, j;
    for (i = 0; i < MAX_STATE; i++)
    {
        initial_sum[i] = 0;
    }
    for (i = 0; i < MAX_STATE; i++)
    {
        for (j = 0; j < MAX_STATE; j++)
        {
            epsilon_sum[i][j] = 0;
            epsilon_sum_div[i][j]=0;
        }
    }
    for (i = 0; i < MAX_SEQ; i++)
    {
        for (j = 0; j < MAX_STATE; j++)
        {
            Gamma_sum[i][j] = 0;
            Gamma_sum_mol[i][j] = 0;
        }
    }
}

void cumulate(HMM *model, char seq[])
{
    int i, j, t, k, seq_len;
    seq_len = strlen(seq);

    // update initial table
    for (i = 0; i < model->state_num; i++)
    {
        initial_sum[i] += Gamma[0][i];
        //model->initial[i] = Gamma[0][i];
    }
    // update transition table
    for (i = 0; i < model->state_num; i++)
    {
        double tmp = 0;
        for (t = 0; t < seq_len - 1; t++)
        {
            tmp += Gamma[t][i];
        }

        for (j = 0; j < model->state_num; j++)
        {
            epsilon_sum[i][j] += epsilon[i][j];
            epsilon_sum_div[i][j]+=tmp;
            //model->transition[i][j] = epsilon[i][j] / tmp;
        }
    }
    // update observation table
    for (i = 0; i < model->state_num; i++)
    {
        for (k = 0; k < model->observ_num; k++)
        {
            double div = 0, mol = 0;
            for (t = 0; t < seq_len; t++)
            {
                div += Gamma[t][i];
                if (seq_to_num(seq[t]) == k)
                {
                    mol += Gamma[t][i];
                }
            }
            Gamma_sum[k][i] += div;
            Gamma_sum_mol[k][i] += mol;
            //model->observation[k][i] = mol / div;
        }
    }
}

void init_epsilon(HMM *model){
    int i,j;
    for (i=0;i<model->state_num;i++){
        for(j=0;j<model->state_num;j++){
            epsilon[i][j]=0;
        }
    }
}

void check_model_valid(HMM *model)
{
    int i, j;
    double initsum = 0;
    for (i = 0; i < model->state_num; i++)
    {
        initsum += model->initial[i];
    }
    printf("init %f\n", initsum);
    for (i = 0; i < model->state_num; i++)
    {
        double sum = 0;
        for (j = 0; j < model->state_num; j++)
        {
            sum += model->transition[i][j];
        }
        printf("A %f\n", sum);
    }
    for (i = 0; i < model->observ_num; i++)
    {
        double sum = 0;
        for (j = 0; j < model->state_num; j++)
        {
            sum += model->observation[j][i];
        }
        printf("B %f\n", sum);
    }
}

void forward(HMM *model, char seq[])
{
    int i, j, t, seq_len;
    seq_len = strlen(seq);
    for (i = 0; i < model->state_num; i++)
    {
        alpha[0][i] = model->observation[seq_to_num(seq[0])][i] * model->initial[i];
    }

    for (t = 1; t < seq_len; t++)
    {
        for (j = 0; j < model->state_num; j++)
        {
            double tmp = 0;
            for (i = 0; i < model->state_num; i++)
            {
                tmp += alpha[t - 1][i] * model->transition[i][j];
            }
            alpha[t][j] = tmp * model->observation[seq_to_num(seq[t])][j];
        }
    }
    /*
    printf("\n This is alpha table:\n");
    for(i=0;i<seq_len;i++){
        for(j=0;j<model->state_num;j++){
            printf("%f ",alpha[i][j]);
        }
        printf("\n");
    }*/
}

void backward(HMM *model, char seq[])
{
    int i, j, t, seq_len;
    seq_len = strlen(seq);
    for (i = 0; i < model->state_num; i++)
    {
        beta[seq_len - 1][i] = 1;
    }
    for (t = seq_len - 2; t >= 0; t--)
    {
        for (i = 0; i < model->state_num; i++)
        {
            beta[t][i] = 0;
            for (j = 0; j < model->state_num; j++)
            {
                beta[t][i] += model->transition[i][j] * model->observation[seq_to_num(seq[t + 1])][j] * beta[t + 1][j];
            }
        }
    } /*
    printf("\n This is beta table:\n");
    for(i=0;i<seq_len;i++){
        for(j=0;j<model->state_num;j++){
            printf("%f ",beta[i][j]);
        }
        printf("\n");
    }*/
}

void BaumWelch(HMM *model, char seq[])
{
    int i, j, t, seq_len;
    seq_len = strlen(seq);
    for (t = 0; t < seq_len; t++)
    {
        for (i = 0; i < model->state_num; i++)
        {
            double div = 0;
            for (j = 0; j < model->state_num; j++)
            {
                div += alpha[t][j] * beta[t][j];
            }
            Gamma[t][i] = (alpha[t][i] * beta[t][i]) / div;
        }
    }

    for (t = 0; t < seq_len - 1; t++)
    {
        double div = 0;
        for (i = 0; i < model->state_num; i++)
        {
            for (j = 0; j < model->state_num; j++)
            {
                div += alpha[t][i] * model->transition[i][j] * model->observation[seq_to_num(seq[t + 1])][j] * beta[t + 1][j];
            }
        }

        for (i = 0; i < model->state_num; i++)
        {
            for (j = 0; j < model->state_num; j++)
            {
                epsilon[i][j] += alpha[t][i] * model->transition[i][j] * model->observation[seq_to_num(seq[t + 1])][j] * beta[t + 1][j] / div;
            }
        }
    } /*
    printf("\n This is gamma table:\n");
    for(i=0;i<seq_len;i++){
        for(j=0;j<model->state_num;j++){
            printf("%f ",Gamma[i][j]);
        }
        printf("\n");
    }
    printf("\n This is epsilon table:\n");
    for(i=0;i<model->state_num;i++){
        for(j=0;j<model->state_num;j++){
            printf("%f ",epsilon[i][j]);
        }
        printf("\n");
    }*/
}

void update_weight(HMM *model, int sen_count, int seq_len)
{
    int i, j, t;
    for (i = 0; i < model->state_num; i++)
    {
        model->initial[i] = initial_sum[i] / sen_count;
    }
    for (i = 0; i < model->state_num; i++)
    {
        for (j = 0; j < model->state_num; j++)
        {
            model->transition[i][j] = epsilon_sum[i][j] / epsilon_sum_div[i][j];
        }
    }
    for (i = 0; i < model->observ_num; i++)
    {
        for (j = 0; j < model->state_num; j++)
        {
            model->observation[i][j] = Gamma_sum_mol[i][j] / Gamma_sum[i][j];
        }
    }
}

void train(HMM *model, char seq_path[], int num_iter)
{
    int iter;
    int first_read = 1;
    char seq[1000];
    FILE *file;
    for (iter = 0; iter < num_iter; iter++)
    {
        int seq_count = 0;
        file = open_or_die(seq_path, "r");
        init_table();
        while (fscanf(file, "%s", seq) != EOF)
        {
            seq_count++;
            init_epsilon(model);
            forward(model, seq);
            backward(model, seq);
            BaumWelch(model, seq);
            cumulate(model, seq);
        }
        update_weight(model, seq_count, strlen(seq));
        fclose(file);
    }
    //check_model_valid(model);
    //dumpHMM(stderr, model);
}

int main(int argc, char *argv[])
{
    // check argv
    if (argc != 5)
    {
        printf("there must be 5 argc\n");
        return 0;
    }

    int num_iter = atoi(argv[1]);
    char *model_init_path = argv[2];
    char *seq_path = argv[3];
    char *output_path = argv[4];

    HMM model;
    loadHMM(&model, model_init_path);
    train(&model, seq_path, num_iter);
    FILE *fileoutput;
    fileoutput = open_or_die(output_path, "w");
    dumpHMM(fileoutput, &model);
    return 0;
}
