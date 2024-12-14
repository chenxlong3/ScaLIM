#ifndef RESULT_INFO_H
#define RESULT_INFO_H
#include "headers.h"
class ResultInfo
{
private:
	double _time_sampling = 0.0, _time_sampling_2 = 0.0,_time_selection=0.0, __Influence = -1.0, __InfluenceOriginal = -1.0, __Approx = -1.0;
	int _k_edges = 0;
	size_t __RRsetsSize = 0;
public:
	ResultInfo()
	{
	}

	~ResultInfo()
	{
	}

	/// Get running time.
	double get_total_time() const
	{
		return this->_time_sampling + this->_time_sampling_2 + _time_selection;
	}
	
	double get_sampling_time() {
		return this->_time_sampling;
	}
	double get_sampling_2_time() {
		return this->_time_sampling_2;
	}

	double get_selection_time() {
		return this->_time_selection;
	}

	/// Get influence spread.
	double get_influence() const
	{
		return __Influence;
	}

	/// Get self-estimated influence spread.
	double get_influence_original() const
	{
		return __InfluenceOriginal;
	}

	/// Get approximation guarantee.
	double get_approximation() const
	{
		return __Approx;
	}

	/// Get seed size.
	int get_k_edges() const
	{
		return _k_edges;
	}

	/// Get the number of RR sets.
	size_t get_RRsets_size() const
	{
		return __RRsetsSize;
	}

	/// Set running time.
	void set_sampling_time(const double value)
	{
		this->_time_sampling = value;
	}
	void set_sampling_2_time(const double value)
	{
		this->_time_sampling_2 = value;
	}
	void set_selection_time(const double value) {
		this->_time_selection = value;
	}
	/// Set influence spread.
	void set_influence(const double value)
	{
		__Influence = value;
	}

	/// Set self-estimated influence spread
	void set_influence_original(const double value)
	{
		__InfluenceOriginal = value;
	}
	
	/// Set approximation guarantee.
	void set_approximation(const double value)
	{
		__Approx = value;
	}

	/// Set seed size.
	void set_k_edges(const int value)
	{
		_k_edges = value;
	}

	/// Set the number of RR sets.
	void set_RR_sets_size(const size_t value)
	{
		__RRsetsSize = value;
	}

	// Save results
	void save_to_file(string filepath) {
		const auto approx = this->get_approximation();
		const auto runTime = this->get_total_time();
		const auto sampling_time = this->get_sampling_time();
		const auto selection_time = this->get_selection_time();
		const auto influence = this->get_influence();
		const auto influenceOriginal = this->get_influence_original();
		
		const auto RRsetsSize = this->get_RRsets_size();
		const auto k_edges = this->get_k_edges();

		std::cout << "   --------------------" << std::endl;
		std::cout << "  |Approx.: " << approx << std::endl;
		std::cout << "  |Time (sec): " << runTime << std::endl;
		std::cout << "  |Sampling Time (sec): " << sampling_time << std::endl;
		std::cout << "  |Selection Time (sec): " << selection_time << std::endl;
		std::cout << "  |Influence: " << influence << std::endl;
		std::cout << "  |Self-estimated influence: " << influenceOriginal << std::endl;
		std::cout << "  |#Edges: " << k_edges << std::endl;
		std::cout << "  |#RR sets: " << RRsetsSize << std::endl;
		std::cout << "   --------------------" << std::endl;
		
		std::ofstream file_stream(filepath);
		if (file_stream.is_open()) {
			file_stream << "Approx.: " << approx << std::endl;
			file_stream << "Time (sec): " << runTime << std::endl;
			file_stream << "Sampling Time (sec): " << sampling_time << std::endl;
			file_stream << "Sampling Phase 2 Time (sec): " << this->_time_sampling_2 << std::endl;
			file_stream << "Selection Time (sec): " << selection_time << std::endl;
			file_stream << "Influence: " << influence << std::endl;
			file_stream << "Self-estimated influence: " << influenceOriginal << std::endl;
			file_stream << "#Edges: " << k_edges << std::endl;
			file_stream << "#RR sets: " << RRsetsSize << std::endl;
			file_stream.close();
		}
    }
};

#endif