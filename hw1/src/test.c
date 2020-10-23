#include "../inc/hmm.h"
#include <stdio.h>

double delta[MAX_SEQ][MAX_STATE];

int seq_to_num(char c)
{
    return c - 'A';
}

double Viterbi(HMM *model, char seq[])
{
    int i, j, t, seq_len;
    seq_len = strlen(seq);
    for (i = 0; i < model->state_num; i++)
    {
        delta[0][i] = model->initial[i] * model->observation[seq_to_num(seq[0])][i];
    }
    for (t = 1; t < seq_len; t++)
    {
        for (j = 0; j < model->state_num; j++)
        {
            double max = 0;
            for (i = 0; i < model->state_num; i++)
            {
                if (delta[t - 1][i] * model->transition[i][j] > max)
                {
                    max = delta[t - 1][i] * model->transition[i][j];
                }
            }
            delta[t][j] = max * model->observation[seq_to_num(seq[t])][j];
        }
    }
    double p = 0;
    for (i = 0; i < model->state_num; i++)
    {
        if (p < delta[seq_len - 1][i])
        {
            p = delta[seq_len - 1][i];
        }
    }
    return p;
}

void test(HMM *models, int model_count, char seq_path[], char output_result_path[])
{
    int i;
    FILE *testfile, *outputfile;
    testfile = open_or_die(seq_path, "r");
    outputfile = open_or_die(output_result_path, "w");
    char seq[1000];
    int flag = 1;
    while (fscanf(testfile, "%s", seq) != EOF)
    {
        if (flag == 0)
        {
            fprintf(outputfile, "\n");
        }
        flag = 0;
        double max_p = 0;
        int max_p_idx;
        for (i = 0; i < model_count; i++)
        {
            double p;
            p = Viterbi(&models[i], seq);
            if (p > max_p)
            {
                max_p = p;
                max_p_idx = i;
            }
        }
        //printf("model %d, probability is %.80f\n", max_p_idx, max_p);
        //printf("model %d\n", max_p_idx);.
        //printf("%s\n", models[max_p_idx].model_name);
        fprintf(outputfile, "%s %e", models[max_p_idx].model_name, max_p);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("There must be 4 argument\n");
    }

    char *models_list_path = argv[1];
    char *seq_path = argv[2];
    char *output_result_path = argv[3];

    int model_count;
    HMM models[10];
    model_count = load_models(models_list_path, models, 10);
    test(models, model_count, seq_path, output_result_path);

    return 0;
}