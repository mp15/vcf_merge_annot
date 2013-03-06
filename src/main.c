// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of vcf merge annot.
//
// vcf merge annot is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#include <htslib/vcf.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>

struct parsed_opts {
    char* input_list_name;
    
    size_t annot_count;
    char** annot_name;

    char* output_name;
};

typedef struct parsed_opts parsed_opts_t;

typedef struct input_list_entry input_list_entry_t;

struct input_list_entry {
    char* name;
    input_list_entry_t* next;
};

struct curr_state {
    input_list_entry_t* input_file_list;
    htsFile* curr_input_file;
    bcf_hdr_t* curr_input_header;

    size_t annot_count;
    htsFile** annot_file;
    bcf1_t** annot_read;
    bcf_hdr_t** annot_header;
    
    htsFile* output_file;
    bcf_hdr_t* output_header;
};

typedef struct curr_state curr_state_t;

parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        printf("Arguments should be: vcf_merge_annot <genotypevcf.list> <sites.vcf> [<sitesX.vcf> ...] <output.vcf>\r\n");
        return NULL;
    }
    
    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;
    
    retval->input_list_name = strdup(argv[1]);

    retval->annot_count = argc-3;
    retval->annot_name = (char**)calloc(retval->annot_count,sizeof(char*));
    size_t i = 0;
    for (; i < retval->annot_count; i++) {
        retval->annot_name[i] = strdup(argv[i+2]);
    }

    retval->output_name = strdup(argv[i+2]);

    return retval;
}

bool load_next_input(curr_state_t* state) {
    if ( state->curr_input_file != NULL ) {
        vcf_close( state->curr_input_file );
    }
    // Move to next entry
    input_list_entry_t* prev = state->input_file_list;
    state->input_file_list = state->input_file_list->next;
    // Clean up
    free(prev->name);
    free(prev);
    // Out of files?
    if ( state->input_file_list == NULL || !strcmp("", state->input_file_list->name) ) return false;
    
    state->curr_input_file = vcf_open(state->input_file_list->name, "r", NULL);
    state->curr_input_header = vcf_hdr_read(state->curr_input_file);
    return true;
}

bool init(parsed_opts_t* opts, curr_state_t** state ) {
    *state = NULL;
    curr_state_t* state_prep = (curr_state_t*) malloc(sizeof(curr_state_t));
    if (state_prep == NULL) return false;
    
    // Open input list
    FILE* input_list = fopen(opts->input_list_name, "r");
    char *line = NULL;
    size_t size = 0;
    state_prep->input_file_list = NULL;
    state_prep->curr_input_file = NULL;
    input_list_entry_t** prev = &(state_prep->input_file_list);
    while (!feof(input_list))
    {
        getline(&line, &size, input_list);
        while (size > 0 && (line[size-1] == '\0' || line[size-1] == '\n')) { line[size-1] = '\0'; size--; } 
        input_list_entry_t* entry = malloc(sizeof(input_list_entry_t));
        entry->name = strndup(line, size);
        entry->next = NULL;
        (*prev) = entry;
        prev = &(entry->next);
    }
    fclose(input_list);
    
    // Load first input
    if ( state_prep->input_file_list == NULL ) return false;
    state_prep->curr_input_file = vcf_open(state_prep->input_file_list->name, "r", NULL);
    state_prep->curr_input_header = vcf_hdr_read(state_prep->curr_input_file);

    // Open files with annotation in
    state_prep->annot_count = opts->annot_count;
    state_prep->annot_file = (vcfFile**)calloc(opts->annot_count, sizeof(vcfFile*));
    state_prep->annot_header = (bcf_hdr_t**)calloc(opts->annot_count, sizeof(bcf_hdr_t*));
    state_prep->annot_read = (bcf1_t**)calloc(opts->annot_count, sizeof(bcf1_t*));
    for (size_t i = 0; i < opts->annot_count; i++) {
        state_prep->annot_file[i] = vcf_open(opts->annot_name[i], "r", NULL);
        if (state_prep->annot_file[i] == NULL) {
            printf("Could not open input file: %s\r\n", opts->annot_name[i]);
            return false;
        }
        state_prep->annot_header[i] = vcf_hdr_read(state_prep->annot_file[i]);
        state_prep->annot_read[i] = bcf_init1();
        vcf_read1(state_prep->annot_file[i], state_prep->annot_header[i], state_prep->annot_read[i]);
        bcf_unpack(state_prep->annot_read[i], BCF_UN_SHR);
    }

    // TODO: consider merging headers instead of just taking first one
    state_prep->output_header = state_prep->curr_input_header;

    state_prep->output_file = vcf_open(opts->output_name, "w", NULL);
    
    if (state_prep->output_file == NULL) {
        printf("Could not open output file: %s\r\n", opts->output_name);
        return false;
    }
    
    *state = state_prep;
    
    return true;
}

bool match(bcf1_t* a, bcf1_t* b)
{
    if (a->rid == b->rid && a->pos == b->pos && a->n_allele == b->n_allele) {
        bool match = true;
        for (int i = 0; i < a->n_allele; i++) {
            if (strcmp(a->d.allele[i], b->d.allele[i])) {
                return false;
            }
        }
        return true;
    }
    return false;
}

// Returns true if position of a greater than b
bool gt(bcf1_t* a, bcf1_t* b)
{
    if (a->rid > b->rid || (a->rid == b->rid && a->pos > b->pos)) {
        return true;
    } else {
        return false;
    }
}

void read_next_annot(curr_state_t* state, int i)
{
    if (vcf_read1(state->annot_file[i], state->annot_header[i], state->annot_read[i]) >= 0)
    {
        bcf_unpack(state->annot_read[i], BCF_UN_SHR);
    }
    else
    {
        bcf_destroy1(state->annot_read[i]);
        state->annot_read[i] = NULL;
        printf("annot file exhausted\n"); // TRACE
    }
}

bool merge(curr_state_t* state) {
    bcf1_t* line = bcf_init1();

    vcf_hdr_write(state->output_file, state->output_header);
    //int output_key = state->output_header.id[BCF_DT_ID].[$key].key;
    do {
        while (vcf_read1(state->curr_input_file, state->curr_input_header, line) >= 0) {
            bcf_unpack(line, BCF_UN_STR);
            for (int i = 0; i < state->annot_count; i++)
            {
                if (state->annot_read[i] != NULL) {
                    if (match(line, state->annot_read[i])) {
                        // copy annots into line
                        printf("match found\n"); //TRACE
                        //line->d.info[output_key].blah = state->annot_read[i]->d.info[state->annot_key[i]].blah;
                        // read next annot
                        read_next_annot(state, i);
                    } else if (gt(line, state->annot_read[i])) { read_next_annot(state, i); }
                }
            }
            vcf_write1(state->output_file, state->output_header, line);
        }
    } while( load_next_input(state) );
    
    bcf_destroy1(line);
    return true;
}

void cleanup(parsed_opts_t* opts, curr_state_t* state) {
    vcf_close(state->output_file);
    free(opts->output_name);
    for (size_t i = 0; i < state->annot_count; i++) {
        vcf_close(state->annot_file[i]);
    }
    for (size_t i = 0; i < opts->annot_count; i++) {
        free(opts->annot_name[i]);
    }

}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    curr_state_t* state = NULL;
    if (!opts || !init(opts, &state)) return -1;
    
    if (!state || !merge(state)) return -1;
    
    cleanup(opts, state);
      
    return 0;
}
