//  * --------------------------------------------------------------------------------------------------------
//  * Name: tspbwt.cpp
//  * Description: A method of Identity by Descent inference from a tree sequence.
//  * Author: Yuan Wei 
//  * Created on: Jan 17, 2025
//  * Modified on: Mar 06, 2025
//  * --------------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <tskit.h>
using namespace std;

//split a string into tokens with given delimiter
vector<string> split_string(string line, char delimiter){
    vector<string> tokens;
    string token;
    stringstream line_stream(line);
    while (getline(line_stream, token, delimiter)){
        if (!token.empty()){
            tokens.push_back(token);
        }
    }
    return tokens;
}

//read genetic map from a file
static map<int, double> read_genetic_map(string input_map_path_and_file_name){
    map<int, double> genetic_map;
    ifstream input_file_data;
    char map_delimiter = '\t';
    input_file_data.open(input_map_path_and_file_name);
    if (!input_file_data){
        cout << "cannot open file " + input_map_path_and_file_name << endl;
        exit(1);
    }
    if (input_file_data.is_open()){
        int line_number = 0;
        string line_from_file;
        while (getline(input_file_data, line_from_file)){
            vector<string> tokens = split_string(line_from_file, map_delimiter);
            if (tokens.size() >= 4){
                //genetic map HapMap format: #Chromosome \t Position(bp) \t Rate(cM/Mb) \t Map(cM)
                if (line_number > 0){
                    //skip the first header line
                    int physical_location = stoi(tokens[1]);
                    double genetic_location = stod(tokens[3]);
                    genetic_map[physical_location] = genetic_location;
                }
            }
            line_number++;
        }
        input_file_data.close();
    }
    return genetic_map;
}

//interpolate the genetic position of a physical position
static double interpolate_genetic_position(map<int, double> & genetic_maps, int physical_position){
    //interpolate genetic position based on physical position from the local tree
    map<int, double>::iterator position_exact_iterator, position_lower_iterator, position_upper_iterator;
    double genetic_position = 0.0;
    //search current physical position from genetic map
    if (!genetic_maps.empty()){
        position_exact_iterator = genetic_maps.find(physical_position);
        if (position_exact_iterator == genetic_maps.end()){
            //exact physical position not found; need to estimate genetic position with neighbor positions in genetic map
            position_upper_iterator = genetic_maps.upper_bound(physical_position);
            if (position_upper_iterator == genetic_maps.end()){
                if (physical_position > genetic_maps.rbegin()->first){
                    //physical position is larger than any physical position in map
                    genetic_position = genetic_maps.rbegin()->second;
                }
                else { 
                    //physical position is less than any physical position in map
                    genetic_position = genetic_maps.begin()->second;
                }
            }
            else {
                if (position_upper_iterator == genetic_maps.begin()){
                    //physical position is less than any physical position in map
                    genetic_position = genetic_maps.begin()->second;
                }
                else {
                    //estimate genetic position by linear interpolation
                    int physical_position_next = position_upper_iterator->first;
                    double genetic_position_next = position_upper_iterator->second;
                    position_lower_iterator = --position_upper_iterator;
                    int physical_position_previous = position_lower_iterator->first;
                    double genetic_position_previous = position_lower_iterator->second;
                    genetic_position = genetic_position_next + (double)(physical_position - physical_position_next) * (genetic_position_next - genetic_position_previous) / (double)(physical_position_next - physical_position_previous);
                }
            }
        }
        else {
            //exact physical position found; just add its corrlated genetic position
            genetic_position = position_exact_iterator->second;
        }
    }
    return genetic_position;
}

//data structure to hold sorted edges of a tree sequence
class SortedEdge{
public:
    double right; //right interval
    tsk_id_t child; //child node id

    SortedEdge(double right_in, tsk_id_t child_in) : right(right_in), child(child_in){
    }

    SortedEdge() : right(0), child(-1){
    }
};

//function to sort edges by their right intervals
bool compare_edges(const SortedEdge &lhs, const SortedEdge &rhs){
    if (lhs.right == rhs.right){
        return lhs.child < rhs.child;
    }
    else {
        return lhs.right < rhs.right;
    }
}

//handle error from tskit library
static void check_tsk_error(int val){
    if (val < 0){
        cerr << "line " << __LINE__ << ": " << tsk_strerror(val) << endl;
        exit(EXIT_FAILURE);
    }
}

//get leaf nodes under a given node
void get_leaf_nodes_under_node(const tsk_tree_t *tree, const tsk_table_collection_t *tables, tsk_id_t curr_node_id, vector<bool> & changed_leaf_nodes, int number_of_samples){
    int ret = 0;
    tsk_size_t num_nodes;
    tsk_id_t *nodes = (tsk_id_t*)tsk_malloc(tsk_tree_get_size_bound(tree) * sizeof(*nodes));
    tsk_id_t v;
    ret = tsk_tree_preorder_from(tree, curr_node_id, nodes, &num_nodes);
    check_tsk_error(ret);
    for (tsk_size_t j = 0; j < num_nodes; j++) {
        v = nodes[j];
        if (v < number_of_samples && tables->nodes.time[v] == 0){
            changed_leaf_nodes[v] = true;
        }
    }
    tsk_safe_free(nodes);
    return;
}

//build pbwt prefix and divergence data
void build_prefix_divergence(int site_index, vector<int> & prefix, vector<int> & divergence, vector<int> & prefix_prev, vector<int> & divergence_prev, vector<bool> & sites){
    int d0t = site_index + 1;
    int d1t = site_index + 1;
    vector<int> prefix_1;
    vector<int> divergence_1;
    for (unsigned long int i = 0; i < sites.size(); i++){
        if (divergence_prev[i] > d0t){
            d0t = divergence_prev[i];
        }
        if (divergence_prev[i] > d1t){
            d1t = divergence_prev[i];
        }
        if (sites[prefix_prev[i]] == 0){
            prefix.emplace_back(prefix_prev[i]);
            divergence.emplace_back(d0t);
            d0t = 0;
        }
        else {
            prefix_1.emplace_back(prefix_prev[i]);
            divergence_1.emplace_back(d1t);
            d1t = 0;
        }
    }
    prefix.insert(prefix.end(), prefix_1.begin(), prefix_1.end());
    divergence.insert(divergence.end(), divergence_1.begin(), divergence_1.end());
    return;
}

//output ibd segments by pbwt for an interval in the middle
void get_long_matches_not_last_site_with_cutoff(int site_index, vector<int> & prefix, vector<int> & divergence, vector<int> & prefix_prev, vector<int> & divergence_prev, vector<bool> & sites, int unit_of_cutoff, double cutoff_length, vector<int> & physical_positions, map<int, double> & genetic_map, unsigned long long int & number_of_ibds, unsigned long long int & length_of_ibds, ofstream & output_data){
    if (cutoff_length > 0){
        vector<int> list_0;
        vector<int> list_1;
        for (unsigned long int i = 0; i < divergence_prev.size(); i++){
            if (unit_of_cutoff == 0 ? physical_positions[site_index + 1] - physical_positions[divergence_prev[i]] < cutoff_length : genetic_map[physical_positions[site_index + 1]] - genetic_map[physical_positions[divergence_prev[i]]] < cutoff_length){
                if (list_0.size() > 0 && list_1.size() > 0){
                    int is = min(list_0[0], list_1[0]);
                    int ie = max(list_0[list_0.size() - 1], list_1[list_1.size() - 1]);
                    for (int i1 = is; i1 <= ie; i1++){
                        int dst = 0;
                        for (int i2 = i1 + 1; i2 <= ie; i2++){
                            if (divergence_prev[i2] > dst){
                                dst = divergence_prev[i2];
                            }
                            if (sites[prefix_prev[i1]] != sites[prefix_prev[i2]]){
                                if (prefix_prev[i1] <= prefix_prev[i2]){
                                    output_data << prefix_prev[i1] << "\t" << prefix_prev[i2] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                                }
                                else {
                                    output_data << prefix_prev[i2] << "\t" << prefix_prev[i1] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                                }
                                number_of_ibds += 1;
                                length_of_ibds += physical_positions[site_index + 1] - physical_positions[dst];
                            }
                        }
                    }
                }
                list_0.clear();
                list_1.clear();
            }
            if (sites[prefix_prev[i]] == 0){
                list_0.emplace_back(i);
            }
            else {
                list_1.emplace_back(i);
            }
        }

        //get matched pairs for the last sample case
        if (list_0.size() > 0 && list_1.size() > 0){
            int is = min(list_0[0], list_1[0]);
            int ie = max(list_0[list_0.size() - 1], list_1[list_1.size() - 1]);
            for (int i1 = is; i1 <= ie; i1++){
                int dst = 0;
                for (int i2 = i1 + 1; i2 <= ie; i2++){
                    if (divergence_prev[i2] > dst){
                        dst = divergence_prev[i2];
                    }
                    if (sites[prefix_prev[i1]] != sites[prefix_prev[i2]]){
                        if (prefix_prev[i1] <= prefix_prev[i2]){
                            output_data << prefix_prev[i1] << "\t" << prefix_prev[i2] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                        }
                        else {
                            output_data << prefix_prev[i2] << "\t" << prefix_prev[i1] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                        }
                        number_of_ibds += 1;
                        length_of_ibds += physical_positions[site_index + 1] - physical_positions[dst];
                    }
                }
            }
        }
        list_0.clear();
        list_1.clear();
    }
    return;
}

//output ibd segments by pbwt for an interval at the end
void get_long_matches_last_site_with_cutoff(int site_index, vector<int> & prefix, vector<int> & divergence, int unit_of_cutoff, double cutoff_length, vector<int> & physical_positions, map<int, double> & genetic_map, unsigned long long int & number_of_ibds, unsigned long long int & length_of_ibds, ofstream & output_data){
    if (cutoff_length > 0){
        vector<int> list_t;
        for (unsigned long int i = 0; i < divergence.size(); i++){
            if (unit_of_cutoff == 0 ? physical_positions[site_index + 1] - physical_positions[divergence[i]] < cutoff_length : genetic_map[physical_positions[site_index + 1]] - genetic_map[physical_positions[divergence[i]]] < cutoff_length){
                if (list_t.size() > 1){
                    for (unsigned long int i1 = 0; i1 < list_t.size(); i1++){
                        int dst = 0;
                        for (unsigned long int i2 = i1 + 1; i2 < list_t.size(); i2++){
                            if (divergence[list_t[i2]] > dst){
                                dst = divergence[list_t[i2]];
                            }
                            if (prefix[list_t[i1]] <= prefix[list_t[i2]]){
                                output_data << prefix[list_t[i1]] << "\t" << prefix[list_t[i2]] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                            }
                            else {
                                output_data << prefix[list_t[i2]] << "\t" << prefix[list_t[i1]] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                            }
                            number_of_ibds += 1;
                            length_of_ibds += physical_positions[site_index + 1] - physical_positions[dst];
                        }
                    }
                }
                list_t.clear();
            }
            list_t.emplace_back(i);
        }

        //report matched sample pair for the last sample case
        if (list_t.size() > 1){
            for (unsigned long int i1 = 0; i1 < list_t.size(); i1++){
                int dst = 0;
                for (unsigned long int i2 = i1 + 1; i2 < list_t.size(); i2++){
                    if (divergence[list_t[i2]] > dst){
                        dst = divergence[list_t[i2]];
                    }
                    if (prefix[list_t[i1]] <= prefix[list_t[i2]]){
                        output_data << prefix[list_t[i1]] << "\t" << prefix[list_t[i2]] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                    }
                    else {
                        output_data << prefix[list_t[i2]] << "\t" << prefix[list_t[i1]] << "\t" << physical_positions[dst] << "\t" << physical_positions[site_index + 1] << endl;
                    }
                    number_of_ibds += 1;
                    length_of_ibds += physical_positions[site_index + 1] - physical_positions[dst];
                }
            }
        }
        list_t.clear();
    }
    return;
}

static void show_usage(string program_name){
    cout << "Usage: " << program_name << " <Option(s)>\n";
    cout << "Option(s):\n";
    cout << "\t-h,--help\t\t\t\t\tShow this help message\n";
    cout << "\t-t,--trees <INPUT TREE SEQUENCE FILE>\t\tInput tree sequence path and file name\n";
    cout << "\t-m,--map <INPUT GENETIC MAP FILE>\t\tInput genetic map path and file name\n";
    cout << "\t-l,--length <MINIMUM CUTOFF LENGTH>\t\tMinimum cutoff length of an IBD (default is 1 base pair)\n";
    cout << "\t-u,--unit <UNIT OF CUTOFF LENGTH>\t\tUnit of cutoff length of an IBD (0: physical; 1: genetic (input genetic map is required); default is 0)\n";
    cout << "\t-o,--output <OUTPUT IBD FILE>\t\t\tOutput ibd path and file name (default is \"./output.ibd\")\n";
}

int main(int argc, char *argv[]){
    try {
        cout << "start program" << endl;

        //variables
        string input_ts_file;
        string input_map_file;
        string output_ibd_file = "./output.ibd";
        double min_cutoff_length = 1; //default is 1 base pair
        int unit_of_cutoff_length = 0; //0: physical; 1: genetic (input genetic map is required); default is 0

        //get command line arguments
        string program_name = argv[0];
        string program_arguments = "";
        for (int i = 1; i < argc; i++){
            string argument = argv[i];
            if ((argument == "-h") || (argument == "--help")){
                show_usage(argv[0]);
                cout << "end program" << endl;
                return 0;
            }
            else if ((argument == "-t") || (argument == "--trees")){
                if (i + 1 < argc){
                    input_ts_file = argv[i + 1];
                    program_arguments += " t=" + input_ts_file;
                    i++;
                }
                else {
                    cout << "-t (or --trees) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-m") || (argument == "--map")){
                if (i + 1 < argc){
                    input_map_file = argv[i + 1];
                    program_arguments += " m=" + input_map_file;
                    i++;
                }
                else {
                    cout << "-m (or --map) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-o") || (argument == "--output")){
                if (i + 1 < argc){
                    output_ibd_file = argv[i + 1];
                    program_arguments += " o=" + output_ibd_file;
                    i++;
                }
                else {
                    cout << "-o (or --output) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-l") || (argument == "--length")){
                if (i + 1 < argc){
                    min_cutoff_length = stod(argv[i + 1]);
                    program_arguments += " l=" + to_string(min_cutoff_length);
                    i++;
                }
                else {
                    cout << "-l (or --length) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else if ((argument == "-u") || (argument == "--unit")){
                if (i + 1 < argc){
                    unit_of_cutoff_length = stoi(argv[i + 1]);
                    program_arguments += " u=" + to_string(unit_of_cutoff_length);
                    i++;
                }
                else {
                    cout << "-u (or --unit) option requires one argument\n";
                    cout << "end program" << endl;
                    return 1;
                }
            }
            else {
                cout << "unrecognized arguments: " << argument << endl;
                show_usage(argv[0]);
                cout << "end program" << endl;
                return 1;
            }
        }
        program_arguments = program_arguments.size() > 0 ? program_arguments.substr(1) : "";

        //output variables
        cout << "program name: " << program_name << "\nprogram arguments: " << program_arguments << endl;

        if (unit_of_cutoff_length != 0 && unit_of_cutoff_length != 1){
            unit_of_cutoff_length = 0;
            cout << "unit_of_cutoff_length is set to: " << unit_of_cutoff_length << endl;
        }

        if (min_cutoff_length <= 0){
            min_cutoff_length = 1;
            cout << "min_cutoff_length is set to: " << min_cutoff_length << endl;
        }

        ofstream output_file_data;
        output_file_data.open(output_ibd_file, ios::trunc);
        if (!output_file_data){
            cout << "cannot create or open file " + output_ibd_file << endl;
            exit(1);
        }
        if (output_file_data.is_open()){
            output_file_data << "#sample_id_1\tsample_id_2\tstart_position\tend_position" << endl;
        }
        
        //load genetic map
        map<int, double> genetic_map;
        if (unit_of_cutoff_length == 1){
            genetic_map = read_genetic_map(input_map_file);
        }

        //load tree sequence
        int ret;
        tsk_treeseq_t ts;
        tsk_table_collection_t tables;
        tsk_tree_t tree;
        ret = tsk_table_collection_load(&tables, input_ts_file.c_str(), 0);
        check_tsk_error(ret);
        ret = tsk_table_collection_build_index(&tables, 0);
        check_tsk_error(ret);
        ret = tsk_treeseq_init(&ts, &tables, 0);
        check_tsk_error(ret);
        ret = tsk_tree_init(&tree, &ts, 0);
        check_tsk_error(ret);
        tsk_size_t num_edges = tsk_treeseq_get_num_edges(&ts);

        //sort edges by right interval
        vector<SortedEdge> sorted_edges;
        sorted_edges.reserve((size_t)num_edges);
        for (size_t e = 0; e < (size_t)num_edges; e++){
            sorted_edges.emplace_back((&tables)->edges.right[e], (&tables)->edges.child[e]);
        }
        sort(begin(sorted_edges), end(sorted_edges), compare_edges);
        
        //initialize variables
        int num_trees = (int)(ts.num_trees);
        int number_of_samples = (int)(ts.num_samples);
        unsigned long int curr_edge_id = 0;
        double start_pos_curr_tree = -1;
        double end_pos_curr_tree = -1;
        map<int, double> ts_genetic_map;
        vector<int> ts_physical_positions;
        unsigned long long int total_number_of_ibds = 0;
        unsigned long long int total_physical_length_of_ibds = 0;
        vector<int> prefix1, prefix2;
        vector<int> divergence1, divergence2;
        vector<int>* prefix_ptr = &prefix1;
        vector<int>* divergence_ptr = &divergence1;
        vector<int>* prefix_prev_ptr = &prefix2;
        vector<int>* divergence_prev_ptr = &divergence2;
        
        //iterate over the tree sequence
        for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)){
            //get tree id and interval positions
            int tree_id = (int)(tree.index);
            start_pos_curr_tree = tree.interval.left;
            end_pos_curr_tree = tree.interval.right;
            
            //update ts physical map
            if (ts_physical_positions.size() == 0 || (ts_physical_positions.size() > 0 && ts_physical_positions[ts_physical_positions.size() - 1] != start_pos_curr_tree)){
                ts_physical_positions.emplace_back(start_pos_curr_tree);
            }
            if (ts_physical_positions.size() > 0 && ts_physical_positions[ts_physical_positions.size() - 1] != end_pos_curr_tree){
                ts_physical_positions.emplace_back(end_pos_curr_tree);
            }

            //update ts genetic map
            if (unit_of_cutoff_length == 1){
                if (ts_genetic_map.count(start_pos_curr_tree) <= 0){
                    ts_genetic_map[start_pos_curr_tree] = interpolate_genetic_position(genetic_map, start_pos_curr_tree);
                }
                if (ts_genetic_map.count(end_pos_curr_tree) <= 0){
                    ts_genetic_map[end_pos_curr_tree] = interpolate_genetic_position(genetic_map, end_pos_curr_tree);
                }
            }

            //get ibd segments
            if (tree_id >= 0 && tree_id < num_trees - 1){
                if (tree.num_edges == 0){
                    //skip the tree
                    if ((*prefix_ptr).size() > 0){
                        get_long_matches_last_site_with_cutoff(tree_id, *prefix_ptr, *divergence_ptr, unit_of_cutoff_length, min_cutoff_length, ts_physical_positions, ts_genetic_map, total_number_of_ibds, total_physical_length_of_ibds, output_file_data);
                    }
                    for (int i = 0; i < number_of_samples; i++){
                        (*prefix_ptr).emplace_back(i);
                        (*divergence_ptr).emplace_back(tree_id + 1);
                    }
                }
                else {
                    //build matrix from sorted edges and run pbwt
                    while (curr_edge_id < sorted_edges.size() && sorted_edges[curr_edge_id].right == end_pos_curr_tree){
                        tsk_id_t c_node_id = sorted_edges[curr_edge_id].child;
                        vector<bool> changed_leaf_nodes(number_of_samples, 0);
                        get_leaf_nodes_under_node(&tree, &tables, c_node_id, changed_leaf_nodes, number_of_samples);
                        if ((*prefix_ptr).size() == 0){
                            for (int i = 0; i < number_of_samples; i++){
                                (*prefix_prev_ptr).emplace_back(i);
                                (*divergence_prev_ptr).emplace_back(0);
                            }
                        }
                        else {
                            if (prefix_prev_ptr == &prefix1){
                                prefix_prev_ptr = &prefix2;
                                divergence_prev_ptr = &divergence2;
                                prefix_ptr = &prefix1;
                                divergence_ptr = &divergence1;
                            }
                            else {
                                prefix_prev_ptr = &prefix1;
                                divergence_prev_ptr = &divergence1;
                                prefix_ptr = &prefix2;
                                divergence_ptr = &divergence2;
                            }
                        }
                        (*prefix_ptr).clear();
                        (*divergence_ptr).clear();
                        build_prefix_divergence(tree_id, *prefix_ptr, *divergence_ptr, *prefix_prev_ptr, *divergence_prev_ptr, changed_leaf_nodes);
                        get_long_matches_not_last_site_with_cutoff(tree_id, *prefix_ptr, *divergence_ptr, *prefix_prev_ptr, *divergence_prev_ptr, changed_leaf_nodes, unit_of_cutoff_length, min_cutoff_length, ts_physical_positions, ts_genetic_map, total_number_of_ibds, total_physical_length_of_ibds, output_file_data);
                        curr_edge_id += 1;
                        changed_leaf_nodes.clear();
                    }
                }
            }
            if (tree_id == num_trees - 1){
                get_long_matches_last_site_with_cutoff(tree_id, *prefix_ptr, *divergence_ptr, unit_of_cutoff_length, min_cutoff_length, ts_physical_positions, ts_genetic_map, total_number_of_ibds, total_physical_length_of_ibds, output_file_data);
            }
        }
        check_tsk_error(ret);

        //clean up
        output_file_data.close();
        sorted_edges.clear();
        (*prefix_ptr).clear();
        (*divergence_ptr).clear();
        (*prefix_prev_ptr).clear();
        (*divergence_prev_ptr).clear();
        tsk_tree_free(&tree);
        tsk_treeseq_free(&ts);
        tsk_table_collection_free(&tables);

        //output summary
        cout << "total_number_of_ibds=" << total_number_of_ibds << endl;
        cout << "total_physical_length_of_ibds=" << total_physical_length_of_ibds << endl;
        cout << "end program" << endl;
    }
    catch (exception & exception_output){
        cerr << exception_output.what() << endl;
    }
    return 0;
}