/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// VQ: Vector quantization

#include <stdlib.h>

#include <cassert>
#include <map>

#include "Basevector.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "system/file/FileReader.h"
#include "system/file/FileWriter.h"

void VQ_normalize(int num_dimensions, vec<String> training_file_list, float* means, float* vars, int RANDOM_SAMPLE_SIZE)
{
    float* sample = new float[num_dimensions];

    RANDOM_SAMPLE_SIZE /= training_file_list.size();
	vec< vec<int> > samples(training_file_list.size());
    for (unsigned int i = 0; i < training_file_list.size(); i++)
	{
        struct stat s;
        stat(training_file_list[i].c_str(), &s);
        int number_of_samples = s.st_size / (sizeof(float) * num_dimensions);

		while (samples[i].size() < (unsigned int)RANDOM_SAMPLE_SIZE)
		{
			int sample_index = rand() % number_of_samples;
			samples[i].push_back(sample_index);
		}
    }

    // First compute the mean of the training set.
    int num_samples = 0;
    float* sum = new float[num_dimensions];
    memset(sum, 0x00, sizeof(float) * num_dimensions);
    for (int i = 0; i < training_file_list.isize(); i++)
    {
        String TRAINING = training_file_list[i];
        printf("Preprocessing (means) %s\n", TRAINING.c_str()); fflush(stdout);
        FileReader training(TRAINING.c_str());

        for (unsigned int sample_index = 0; sample_index < samples[i].size(); sample_index++)
        {
            training.seek(sizeof(float)*num_dimensions*samples[i][sample_index]);
            training.read(sample,num_dimensions*sizeof(float));
            for (int j = 0; j < num_dimensions; j++)
            {
                sum[j] += sample[j];
            }
            num_samples += 1;
        }
    }
    for (int j = 0; j < num_dimensions; j++)
    {
        means[j] = sum[j] / (float)num_samples;
    }

    // Now compute the standard deviations.
    memset(sum, 0x00, sizeof(float) * num_dimensions);
    num_samples = 0;
    for (int i = 0; i < training_file_list.isize(); i++)
    {
        String TRAINING = training_file_list[i];
        printf("Preprocessing (std) %s\n", TRAINING.c_str()); fflush(stdout);
        FileReader training(TRAINING.c_str());

        for (unsigned int sample_index = 0; sample_index < samples[i].size(); sample_index++)
        {
            training.seek(sizeof(float)*num_dimensions*samples[i][sample_index]);
            training.read(sample,num_dimensions*sizeof(float));

            for (int j = 0; j < num_dimensions; j++)
            {
                sum[j] += (sample[j] - means[j])*(sample[j] - means[j]);
            }
            num_samples += 1;
        }
    }
    for (int j = 0; j < num_dimensions; j++)
    {
        vars[j] = sqrt(sum[j] / (float)num_samples);
        if (vars[j] == 0.0) { vars[j] = 1.00000000; } // This must be a binary feature.
    }

    // The class labels are always at index 0.
    means[0] = 0;
    vars[0]  = 1;
}

// Randomly choose codebook vectors.
void VQ_init_standard(int codebook_size, float*** codebook, int num_dimensions, vec<String> training_file_list, float* means, float* vars)
{
    *codebook = new float*[codebook_size];
    float* sample = new float[num_dimensions];

    // Now randomly select data.
    for (int i  = 0; i < codebook_size; i++)
    {
        (*codebook)[i] = new float[num_dimensions+1];

        // Randomly pick a sample from the training set.
        int training_file_index = rand() % training_file_list.size();
        struct stat s;
        stat(training_file_list[training_file_index].c_str(), &s);
        int training_file_size = s.st_size / (sizeof(float) * num_dimensions);
        int training_file_record = rand() % training_file_size;
        FileReader f(training_file_list[training_file_index].c_str() );
        f.seek(training_file_record * (sizeof(float)*num_dimensions));
        f.read(sample, num_dimensions*sizeof(float));

        float radius = 0;
        for (int j = 0; j < num_dimensions; j++)
        {
            (*codebook)[i][j] = (sample[j] - means[j]) / vars[j];
            radius += (*codebook)[i][j] * (*codebook)[i][j];
        }
        radius = sqrt(radius);
        printf("codebook radius: %f\n", radius);
    }
}

void VQ_init_balanced_radius(int codebook_size, float*** codebook, int num_dimensions, vec<String> training_file_list, float* means, float* vars)
{
    *codebook = new float*[codebook_size];
    float* sample = new float[num_dimensions];

    int sample_set_1_size = 10000;
    vec< pair<int,int> > sample_set_1; sample_set_1.reserve(sample_set_1_size);
    vec<float> radii; 

    // Randomly select a LOT of samples.
    for (int i  = 0; i < sample_set_1_size; i++)
    {
        int training_file_index = rand() % training_file_list.size();
        struct stat s;
        stat(training_file_list[training_file_index].c_str(), &s);
        int training_file_size = s.st_size / (sizeof(float) * num_dimensions);
        int training_file_record = rand() % training_file_size;
        FileReader f(training_file_list[training_file_index].c_str());
        f.seek(training_file_record * (sizeof(float)*num_dimensions));
        f.read(sample, sizeof(float)*num_dimensions);

        float radius = 0.0;
        for (int j = 1; j < num_dimensions; j++)
        {
            sample[j] = (sample[j] - means[j]) / vars[j];
            radius += sample[j] * sample[j];
        }
        radius = sqrt(radius);

        pair<int,int> p(training_file_index, training_file_record);
        sample_set_1.push_back(p);
        radii.push_back(radius);
    }

    // Draw a population of samples evenly distributed across the radius range.
    vec<float> radii_sorted = radii;
    vec< pair<int,int> > sample_set_1_sorted = sample_set_1;
    SortSync(radii_sorted, sample_set_1_sorted);

    vec< pair<int,int> > sample_set_2;
    int divisions = 4;
    for (int quartile = 0; quartile < divisions; quartile += 1)
    {
        for (int i = 0; i < codebook_size/divisions; i++)
        {
            int s = (quartile*(radii_sorted.size()/divisions)) + (rand() % (radii_sorted.size()/divisions));
            sample_set_2.push_back(sample_set_1_sorted[s]);
        }
    }

    assert(sample_set_2.isize() == codebook_size);

    // Initialize the codebook.
    for (int i  = 0; i < codebook_size; i++)
    {
        int training_file_index  = sample_set_2[i].first;
        int training_file_record = sample_set_2[i].second;

        FileReader f(training_file_list[training_file_index].c_str());
        f.seek(training_file_record * (sizeof(float)*num_dimensions));
        f.read(sample, sizeof(float)*num_dimensions);

        (*codebook)[i] = new float[num_dimensions+1];

        float radius = 0;
        for (int j = 1; j < num_dimensions; j++)
        {
            (*codebook)[i][j] = (sample[j] - means[j]) / vars[j];
            radius += (*codebook)[i][j] * (*codebook)[i][j];
        }
        radius = sqrt(radius);
        printf("codebook radius: %f\n", radius);
    }
    
}

void VQ_load(const char* file_name, int* codebook_size, float*** codebook, int sample_size, float** labels)
{
    printf("Loading %s\n", file_name);

    FileReader codebook_file(file_name);
    float sample[512];
    float y; 
    memset(sample, 0x00, sizeof(float)*512);

    *codebook_size = 0;
    size_t len = sample_size*sizeof(float);
    while ( codebook_file.readSome(sample,len) == len )
    {
        (*codebook_size)++; 
    }

    *codebook = new float*[*codebook_size];
    *labels = new float[*codebook_size];

    codebook_file.seek(0);

    int i = 0;
    while (codebook_file.readSome(sample,len) == len )
    {
        (*codebook)[i] = new float[sample_size];
        memcpy((*codebook)[i], sample, sizeof(float) * (sample_size));
        (*labels)[i] = (*codebook)[i][0];
        (*codebook)[i][0] = 0.0f;
        i++;
    }
}

int VQ_nearest(float** codebook, int codebook_size, float* sample, int sample_size)
{
    int    min_index    = -1;
    double min_distance = 1000000000000.0;

    for (int i = 0; i < codebook_size; i++)
    {
        double distance = 0.0;
        for (int j = 1; j < sample_size; j++)
        {
            distance += pow(codebook[i][j] - sample[j],2);
        }
        distance = sqrt(distance);
#if 0
        printf("CODEBOOK %d:", i);
        for (int j = 0; j < sample_size; j++)
        {
            printf("%f ", codebook[i][j]);
        }
        printf("\n"); fflush(stdout);

        printf("SAMPLE: ");
        for (int j = 0; j < sample_size; j++)
        {
            printf("%f ", sample[j]);
        }
        printf("\n"); fflush(stdout);

        printf(">> %d %f %f\n", i, distance, min_distance);
#endif

        if (distance < min_distance)
        {
            min_distance = distance;
            min_index    = i;
        }
    } 

    assert(min_index != -1);
    
    return min_index;
}

int VQ_update(float** codebook, int codebook_size, float* sample, int sample_size, double learning_rate)
{
    int min_index;
    min_index = VQ_nearest(codebook, codebook_size, sample, sample_size);

    // Update the winning codebook.
    if (learning_rate != 0.0)
    {
        for (int i = 1; i < sample_size; i++)
        {
            double delta = (sample[i] - codebook[min_index][i]) * learning_rate;
            codebook[min_index][i] += delta;
        } 
    }

    return min_index;
}


int VQ_training(vec<String>          training_file_list, 
                vec<String>          calibration_file_list, 
                String               MODEL, 
                int                  NUM_FEATURES, 
                int                  NUM_CODEBOOKS, 
                float                LEARNING_RATE, 
                float                COOLING_RATE, 
                int                  NUM_ITERATIONS, 
                float                SUPERVISION_WEIGHT, 
                String               MEANS_AND_VARS, 
                bool                 LOAD_MEANS_AND_VARS, 
                int                  RANDOM_SAMPLE_SIZE, 
                bool                 BALANCED_RADIUS_INIT, 
                int                  MAX_DEPTH, 
                int                  MIN_SAMPLE_SIZE, 
                String               VQ_EXECUTABLE_NAME)
{
    float** codebook = NULL;
    int sample_size = NUM_FEATURES;
    float sample[512];
    int y;

    vec<int> num_error(NUM_CODEBOOKS);
    vec<int> counts(NUM_CODEBOOKS);

    float* means = new float[NUM_FEATURES];
    float* vars  = new float[NUM_FEATURES];

    if (LOAD_MEANS_AND_VARS)
    {
        FileReader means_and_vars_file(MEANS_AND_VARS.c_str());
        means_and_vars_file.read(means, sizeof(float)*NUM_FEATURES);
        means_and_vars_file.read(vars, sizeof(float)*NUM_FEATURES);
    }
    else
    {
	    VQ_normalize(sample_size, training_file_list, means, vars, RANDOM_SAMPLE_SIZE);
	    FileWriter means_and_vars_file(MEANS_AND_VARS.c_str());
	    means_and_vars_file.write(means, sizeof(float)*NUM_FEATURES);
	    means_and_vars_file.write(vars,  sizeof(float)*NUM_FEATURES);
    }

    if (BALANCED_RADIUS_INIT)
    {
        VQ_init_balanced_radius(NUM_CODEBOOKS, &codebook, sample_size, training_file_list, means, vars);
    }
    else
    {
        VQ_init_standard(NUM_CODEBOOKS, &codebook, sample_size, training_file_list, means, vars);
    }

   // Find the bins.
   for (int training_index = 0; training_index < training_file_list.isize(); training_index += 1)
    {
        String TRAINING = training_file_list[training_index];

        printf("Training from %s\n", TRAINING.c_str()); fflush(stdout);

        for (int t = 0; t < NUM_ITERATIONS; t++)
        {
            FileReader training(TRAINING.c_str());

            if (LEARNING_RATE <= 0.000001) { break; }

            printf("ITERATION %d LEARNING_RATE %f COOLING_RATE %f\n",   
                    t, LEARNING_RATE, COOLING_RATE); fflush(stdout);

            int sampleindex = 0;
            size_t len = sample_size*sizeof(float);
            while ( training.readSome(sample,len) == len )
            { 
                for (int j = 0; j < sample_size; j++)
                {
                    sample[j] = (sample[j] - means[j]) / vars[j];
                }

                y = (int)sample[0];

                sample[0] *= SUPERVISION_WEIGHT;

                if (sampleindex % (1024*1024) == 0)
                {
                    printf("(PASS 1) SAMPLE %d               \r", sampleindex); fflush(stdout);
                }
                sampleindex += 1;

                VQ_update(codebook, NUM_CODEBOOKS, sample, sample_size, LEARNING_RATE);
            }

            LEARNING_RATE *= COOLING_RATE;
        }
    }


        // Initalize the output sample lists.
        vec< vec< pair<int,int> > > sample_lists(NUM_CODEBOOKS);
        for (int i = 0; i < NUM_CODEBOOKS; i++)
        {
            sample_lists[i].reserve(20000000);
            sample_lists[i].resize(0);
        }


    // Estimate the quality for each bin.
   for (int calibration_index = 0; calibration_index < calibration_file_list.isize(); calibration_index += 1)
    {
        String CALIBRATION = calibration_file_list[calibration_index];

        printf("Calibrating from %s\n", CALIBRATION.c_str()); fflush(stdout);

        FileReader calibration(CALIBRATION.c_str());
        int sampleindex = 0;

        size_t len = sample_size*sizeof(float);
        while ( calibration.readSome(sample,len) == len )
        {
                for (int j = 0; j < sample_size; j++)
                {
                    sample[j] = (sample[j] - means[j]) / vars[j];
                }

                y = (int)sample[0];
                sample[0] *= SUPERVISION_WEIGHT;

                if (sampleindex % (1024*1024) == 0)
                {
                    printf("(PASS 2) SAMPLE %d               \r", sampleindex); fflush(stdout);
                }
                sampleindex += 1;

            int nearest_codebook = VQ_nearest(codebook, NUM_CODEBOOKS, sample, sample_size);

                    pair<int,int> sample_id(calibration_index, sampleindex);
                    sample_lists[nearest_codebook].push_back(sample_id);

            if (y == -1) { num_error[nearest_codebook] += 1; }
            counts[nearest_codebook] += 1;
        }
    }

   // Output the model
   FileWriter model(MODEL.c_str());
   for (int i = 0; i < NUM_CODEBOOKS; i++)
   {
	   if (counts[i] == 0) { continue; }
	   float Q = -10.0 * log10((double)num_error[i] / (double)counts[i]);

       codebook[i][0] = Q;
       model.write(codebook[i], sizeof(float)*sample_size);
   }
   model.close();
    
   // Output the model bin sizes
   String MODEL_BIN_SIZES = MODEL + ".bin_sizes";
   std::ofstream model_bin_sizes(MODEL_BIN_SIZES.c_str());
   model_bin_sizes << std::fixed << std::setprecision(1);
   for (int i = 0; i < NUM_CODEBOOKS; i++)
   {
	   if (counts[i] == 0) { continue; }
	   double Q = -10 * log10((double)num_error[i] / (double)counts[i]);
	   model_bin_sizes << Q << " ( " << num_error[1] << " / "
	                   << counts[i] << " ) \n";
   }
   model_bin_sizes.close();
   ForceAssert(model_bin_sizes);

    // If we're doing a hierarchical model...
    if (MAX_DEPTH > 1)
    {
        // for each codebook...
        for (int i = 0; i < NUM_CODEBOOKS; i++)
        {
            // ... if our sample size is too small, don't recurse.
            if (counts[i] < MIN_SAMPLE_SIZE) { continue; }

            // ... otherwise recurse!

            // figure out the parameters,
            String sample_file_name         = MODEL + "." + ToString(i) + ".features";
            String model_file_name          = MODEL + "." + ToString(i);
            String means_and_vars_file_name = MEANS_AND_VARS + "." + ToString(i);

            String training_command = VQ_EXECUTABLE_NAME + 
                                        " MODEL="                + model_file_name                    + 
                                        " TRAINING="             + sample_file_name                   +
                                        " CALIBRATION="          + sample_file_name                   +
                                        " NUM_FEATURES="         + ToString(NUM_FEATURES)             + 
                                        " NUM_CODEBOOKS="        + ToString(NUM_CODEBOOKS)            + 
                                        " LEARNING_RATE="        + ToString(LEARNING_RATE)            + 
                                        " COOLING_RATE="         + ToString(COOLING_RATE)             + 
                                        " NUM_ITERATIONS="       + ToString(NUM_ITERATIONS)           + 
                                        " SUPERVISION_WEIGHT="   + ToString(SUPERVISION_WEIGHT)       + 
                                        " MEANS_AND_VARS="       + means_and_vars_file_name           +
                                        " LOAD_MEANS_AND_VARS="  + (LOAD_MEANS_AND_VARS ? "True" : "False") +
                                        " RANDOM_SAMPLE_SIZE="   + ToString(RANDOM_SAMPLE_SIZE)       + 
                                        " BALANCED_RADIUS_INIT=" + (BALANCED_RADIUS_INIT ? "True" : "False")     + 
                                        " MAX_DEPTH="            + ToString(MAX_DEPTH-1)              + //Note that we drop the depth with each recursion.
                                        " MIN_SAMPLE_SIZE="      + ToString(MIN_SAMPLE_SIZE)          ; 

            // write out the sample file,
            printf("Writing %s\n", sample_file_name.c_str()); fflush(stdout);

            vec<FileReader*> fds(training_file_list.size());
            for (int idx = 0; idx < training_file_list.isize(); idx++)
                fds[idx] = new FileReader(training_file_list[idx].c_str());

            FileWriter sample_file(sample_file_name.c_str());
            for (int j = 0; j < sample_lists[i].isize(); j++)
            {
                int training_file_index = sample_lists[i][j].first;
                FileReader* pFR = fds[training_file_index];
                int sample_index = sample_lists[i][j].second;

                pFR->seek(sample_index*sizeof(float)*NUM_FEATURES);
                pFR->read(sample, sizeof(float)*NUM_FEATURES);
                sample_file.write(sample, sizeof(float)*NUM_FEATURES);
            }

            for (int idx = 0; idx < training_file_list.isize(); idx++)
                delete fds[idx];

            sample_file.close();
 
            // then run the trainer.
            SystemSucceed(training_command);

            // Then delete the features to keep disk usage bounded.
            unlink(sample_file_name.c_str()); 
        }
    }


   return 1;
}

// Manage a folder full of hierarchical VQ models for classification 
// -- because we can't just recurse during testing.
class VQ_model_manager
{

public:

    VQ_model_manager(String MODEL, String MEANS_AND_VARS, int NUM_FEATURES, int MIN_SAMPLE_SIZE)
    {
        this->num_features    = NUM_FEATURES;
        this->top_model       = MODEL;
        this->min_sample_size = MIN_SAMPLE_SIZE;
        this->add(MODEL, MEANS_AND_VARS, NUM_FEATURES);
    }

    // Recursively load each a codebook and all of its children.
    bool add(String model, String means_and_vars, int NUM_FEATURES)
    {
        int codebook_size;
        float** codebook;
        float*  labels;
        float*  means;
        float*  vars;

        VQ_load(model.c_str(), &codebook_size, &codebook, NUM_FEATURES, &labels);
        means = new float[NUM_FEATURES];
        vars  = new float[NUM_FEATURES];
        if ( true )
        {
            FileReader stats(means_and_vars.c_str());
            stats.read(means, sizeof(float)*NUM_FEATURES);
            stats.read(vars,  sizeof(float)*NUM_FEATURES);
        }

        this->codebooks[model] = codebook;
        this->sizes[model]     = codebook_size;
        this->labels[model]    = labels;
        this->means[model]     = means;
        this->vars[model]      = vars;

        Ifstream(bin_sizes_file, model + ".bin_sizes");
        vec<int> bin_sizes(codebook_size);
        for (int i = 0; i < codebook_size; i++)
        {
            String s;
            getline(bin_sizes_file, s);
            vec<String> tokens;
            Tokenize(s, tokens);
            int bin_size = tokens[4].Int();
            bin_sizes[i] = bin_size;
        }
        this->bin_sizes[model] = bin_sizes;

        for (int i = 0; i < codebook_size; i++)
        {
            String sub_model          = model          + "." + ToString(i);
            String sub_means_and_vars = means_and_vars + "." + ToString(i);
            if (IsRegularFile(sub_model)) 
            {
                this->add(sub_model, sub_means_and_vars, NUM_FEATURES);
            }
        }

        return true;
    }
    
    float classify(float* sample)
    {
        float q;
        this->classify(top_model, sample, &q);
        return q;
    }

    bool classify(String model, float* sample, float* q)
    {
        float* means = this->means[model];
        float* vars  = this->vars[model];
        float  local_sample[512];

        //printf("CLASSIFY: %d %s\n", this->num_features, model.c_str()); fflush(stdout);
        for (int i = 0; i < this->num_features; i++)
        {
            local_sample[i] = (sample[i] - means[i]) / vars[i];
            //printf("SAMPLE %i %f %f %f %f\n", i, local_sample[i], sample[i], means[i], vars[i]);
        }

        int   local_v = VQ_nearest(codebooks[model], sizes[model], local_sample, this->num_features);
        float local_q = labels[model][local_v];

        //printf("LOCAL %d %f\n", local_v, local_q); fflush(stdout);

        // If we're below the minimum sample size, return false.
        if (bin_sizes[model][local_v] < this->min_sample_size) { return false; }

        // If we're not at the bottom, recurse.
        String sub_model = model + "." + ToString(local_v);
        if (codebooks.count(sub_model) > 0)
        {
            //printf("RECURSING.\n"); fflush(stdout);
            bool status = classify(sub_model, sample, q);
            if (status) 
            {
                *q = Max(*q, local_q);
                return true;
            }
            else
            {
                *q = local_q;
                return true;
            }
        }
        else
        {
            *q = local_q;
            return true;
        }
    }

private:

    String top_model;
    int num_features;
    int min_sample_size;

    map<String, float**>   codebooks;
    map<String, int>       sizes;
    map<String, vec<int> > bin_sizes;
    map<String, float*>    labels;
    map<String, float*>    means;
    map<String, float*>    vars;

};

int VQ_test(vec<String> test_file_list, String MODEL, int NUM_FEATURES, int NUM_CODEBOOKS, float LEARNING_RATE, float COOLING_RATE, int NUM_ITERATIONS, float SUPERVISION_WEIGHT, String QUALITIES_OUTPUT, String MEANS_AND_VARS, int MIN_SAMPLE_SIZE)
{
        char buf[8192];

        VQ_model_manager model(MODEL, MEANS_AND_VARS, NUM_FEATURES, MIN_SAMPLE_SIZE);

        for (int test_file_index = 0; test_file_index < test_file_list.isize(); test_file_index += 1)
        {
	        String TEST = test_file_list[test_file_index];

            printf("TEST FILE : %s\n", TEST.c_str()); fflush(stdout);
	
	        FileReader test_file(TEST.c_str());
	        std::ofstream out(QUALITIES_OUTPUT.c_str());
	        int sample_size = NUM_FEATURES;
	        float sample[512];
	        int y;
	        int sample_number = 0;
	        printf("sample_size: %d features.\n\n", sample_size); fflush(stdout);
	        size_t len = sizeof(float)*sample_size;
		while ( test_file.readSome(sample,len) == len )
	        {
	            y = (int)sample[0];
	            sample[0] *= SUPERVISION_WEIGHT;
	
	            float quality = model.classify(sample);
	            out << quality << " NA " << y << '\n';
	
	            if (sample_number % (1024*1024) == 0) { printf("processed %d bases.         \r", sample_number); fflush(stdout); }
	            sample_number++;
	        }
		out.close();
		ForceAssert(out);
        }
        
        return 1;
}

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(MODEL);

  CommandArgument_Int_OrDefault(NUM_FEATURES, 11);

  CommandArgument_String_OrDefault(TRAINING, "");
  CommandArgument_String_OrDefault(CALIBRATION, "");

  CommandArgument_String_OrDefault(TRAINING_LIST, "");
  CommandArgument_String_OrDefault(CALIBRATION_LIST, "");

  CommandArgument_Int_OrDefault(NUM_CODEBOOKS,    20);
  CommandArgument_Double_OrDefault(LEARNING_RATE, 0.3);
  CommandArgument_Double_OrDefault(COOLING_RATE,  0.1);
  CommandArgument_Int_OrDefault(NUM_ITERATIONS,   1);
  CommandArgument_Double_OrDefault(SUPERVISION_WEIGHT,  0.0);

  CommandArgument_String_OrDefault(TEST_LIST, "");
  CommandArgument_String_OrDefault(TEST, "");

  CommandArgument_String_OrDefault(QUALITIES_OUTPUT, "");

  CommandArgument_String_OrDefault(MEANS_AND_VARS, "");
  CommandArgument_Bool_OrDefault(LOAD_MEANS_AND_VARS, false);
  CommandArgument_Int_OrDefault(RANDOM_SAMPLE_SIZE, 100000);
  CommandArgument_Bool_OrDefault(BALANCED_RADIUS_INIT, false);

    // parameters for hierarchical VQ
    CommandArgument_Int_OrDefault(MAX_DEPTH, 1);
    CommandArgument_Int_OrDefault(MIN_SAMPLE_SIZE, 1000000);
    CommandArgument_Int_OrDefault(MIN_SAMPLE_SIZE_TEST, 1000000);

  EndCommandArguments;

    String VQ_EXECUTABLE_NAME = argv[0];

    if ((TRAINING != "") || (TRAINING_LIST != ""))
    {
        vec<String> training_file_list;
        vec<String> calibration_file_list;

        if (TRAINING_LIST != "")
        {
            Ifstream(training_list, TRAINING_LIST);
            while (! training_list.eof() )
            {
                String s;
                training_list >> s;
                if (training_list.eof()) { break; }
                training_file_list.push_back(s);
            }
            training_list.close();
        }
        else
        {
            training_file_list.push_back(TRAINING);
        }

        if (CALIBRATION_LIST != "")
        {
            Ifstream(calibration_list, CALIBRATION_LIST);
            while (! calibration_list.eof() )
            {
                String s;
                calibration_list >> s;
                if (calibration_list.eof()) { break; }
                calibration_file_list.push_back(s);
            }
            calibration_list.close();
        }
        else
        {
            calibration_file_list.push_back(CALIBRATION);
        }

        VQ_training(training_file_list, 
                    calibration_file_list, 
                    MODEL, 
                    NUM_FEATURES, 
                    NUM_CODEBOOKS, 
                    LEARNING_RATE, 
                    COOLING_RATE, 
                    NUM_ITERATIONS, 
                    SUPERVISION_WEIGHT, 
                    MEANS_AND_VARS, 
                    LOAD_MEANS_AND_VARS, 
                    RANDOM_SAMPLE_SIZE, 
                    BALANCED_RADIUS_INIT,
                    MAX_DEPTH,
                    MIN_SAMPLE_SIZE,
                    VQ_EXECUTABLE_NAME);
    }

    if ((TEST != "") || (TEST_LIST != ""))
    {
        vec<String> test_file_list;

        if (TEST_LIST != "")
        {
            Ifstream(test_list, TEST_LIST);
            while (! test_list.eof() )
            {
                String s;
                test_list >> s;
                if (test_list.eof()) { break; }
                test_file_list.push_back(s);
            }
            test_list.close();
        }
        else
        {
            test_file_list.push_back(TEST);
        }

        VQ_test(test_file_list, MODEL, NUM_FEATURES, NUM_CODEBOOKS, LEARNING_RATE, COOLING_RATE, NUM_ITERATIONS, SUPERVISION_WEIGHT, QUALITIES_OUTPUT, MEANS_AND_VARS, MIN_SAMPLE_SIZE_TEST);
    }
}
