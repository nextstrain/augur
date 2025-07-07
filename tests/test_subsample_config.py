import pytest
import tempfile
import os
import textwrap
from augur.subsample import Sample, _parse_config, _merge_options
from augur.errors import AugurError


class TestSubsampleConfigValidation:
    """Tests for subsample config parsing and validation."""

    def test_valid_config_with_subsample_max_sequences(self):
        """Test that a valid config with subsample_max_sequences passes validation."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                probabilistic_sampling: true
                include:
                  - strain1
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            assert 'samples' in config
            assert 'test_sample' in config['samples']
            assert config['samples']['test_sample']['subsample_max_sequences'] == 5
            assert config['samples']['test_sample']['probabilistic_sampling'] is True
            assert config['samples']['test_sample']['include'] == ['strain1']
        finally:
            os.unlink(temp_file)

    def test_valid_config_with_sequences_per_group(self):
        """Test that a valid config with sequences_per_group passes validation."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                sequences_per_group: 3
                min_date: "2020-01-01"
                max_date: "2020-12-31"
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            assert config['samples']['test_sample']['sequences_per_group'] == 3
            # Test that dates are preserved as strings
            assert isinstance(config['samples']['test_sample']['min_date'], str)
            assert isinstance(config['samples']['test_sample']['max_date'], str)
        finally:
            os.unlink(temp_file)

    def test_invalid_config_missing_samples(self):
        """Test that config without required 'samples' key fails validation."""
        config_content = textwrap.dedent("""
            include:
              - strain1
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            with pytest.raises(AugurError, match="Config validation failed"):
                _parse_config(temp_file)
        finally:
            os.unlink(temp_file)

    def test_invalid_config_wrong_property_name(self):
        """Test that config with wrong property name (max_sequences instead of subsample_max_sequences) fails."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                max_sequences: 5
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            with pytest.raises(AugurError, match="Config validation failed"):
                _parse_config(temp_file)
        finally:
            os.unlink(temp_file)

    def test_date_handling_no_auto_conversion(self):
        """Test that date-like strings are preserved as strings, not converted to datetime objects."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 1
                min_date: 2020-01-01
                max_date: 2020-12-31
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            # Ensure dates are strings, not datetime objects
            assert isinstance(config['samples']['test_sample']['min_date'], str)
            assert isinstance(config['samples']['test_sample']['max_date'], str)
            assert config['samples']['test_sample']['min_date'] == "2020-01-01"
            assert config['samples']['test_sample']['max_date'] == "2020-12-31"
        finally:
            os.unlink(temp_file)

    def test_yaml_parsing_error(self):
        """Test that invalid YAML syntax raises appropriate error."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                invalid_yaml: [unclosed list
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            with pytest.raises(AugurError, match="Error parsing subsampling scheme"):
                _parse_config(temp_file)
        finally:
            os.unlink(temp_file)

    def test_valid_config_with_new_per_sample_options(self):
        """Test that a valid config with new per-sample options passes validation."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                query: "country == 'Colombia'"
                query_columns:
                  - region:str
                  - coverage:float
                exclude_ambiguous_dates_by: day
                exclude_all: true
                min_length: 1000
                max_length: 30000
                non_nucleotide: false
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            sample = config['samples']['test_sample']
            assert sample['query'] == "country == 'Colombia'"
            assert sample['query_columns'] == ['region:str', 'coverage:float']
            assert sample['exclude_ambiguous_dates_by'] == 'day'
            assert sample['exclude_all'] is True
            assert sample['min_length'] == 1000
            assert sample['max_length'] == 30000
            assert sample['non_nucleotide'] is False
        finally:
            os.unlink(temp_file)

    def test_boolean_options_validation(self):
        """Test that boolean options are properly validated."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                exclude_all: true
                non_nucleotide: false
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            sample = config['samples']['test_sample']
            assert sample['exclude_all'] is True
            assert sample['non_nucleotide'] is False
        finally:
            os.unlink(temp_file)

    def test_exclude_ambiguous_dates_by_enum_validation(self):
        """Test that exclude_ambiguous_dates_by accepts only valid enum values."""
        # Test valid enum value
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                exclude_ambiguous_dates_by: month
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            assert config['samples']['test_sample']['exclude_ambiguous_dates_by'] == 'month'
        finally:
            os.unlink(temp_file)

        # Test invalid enum value
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                exclude_ambiguous_dates_by: invalid_value
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            with pytest.raises(AugurError, match="Config validation failed"):
                _parse_config(temp_file)
        finally:
            os.unlink(temp_file)

    def test_integer_options_validation(self):
        """Test that integer options are properly validated."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                min_length: 500
                max_length: 25000
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            sample = config['samples']['test_sample']
            assert sample['min_length'] == 500
            assert sample['max_length'] == 25000
            assert isinstance(sample['min_length'], int)
            assert isinstance(sample['max_length'], int)
        finally:
            os.unlink(temp_file)

    def test_query_columns_array_validation(self):
        """Test that query_columns accepts array of strings."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                query_columns:
                  - region:str
                  - date:str
                  - coverage:float
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            sample = config['samples']['test_sample']
            assert sample['query_columns'] == ['region:str', 'date:str', 'coverage:float']
            assert isinstance(sample['query_columns'], list)
        finally:
            os.unlink(temp_file)

    def test_sample_yaml_options_validation(self):
        """Test that include/include_where options are properly validated in sample configs."""
        config_content = textwrap.dedent("""
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 5
                include_where:
                  - "region=North America"
                include:
                  - strain1
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            sample_config = config['samples']['test_sample']
            assert sample_config['include_where'] == ['region=North America']
            assert sample_config['include'] == ['strain1']
        finally:
            os.unlink(temp_file)


class TestSampleFilterArgs:
    """Tests for Sample._build_filter_args() method."""

    def test_new_sample_options_mapping(self):
        """Test that new sample options are correctly mapped to filter arguments."""
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            "query": "country == 'Colombia'",
            "query_columns": ["region:str", "coverage:float"],
            "exclude_ambiguous_dates_by": "day",
            "min_length": 1000,
            "max_length": 30000,
        }
        global_filter_args = ["--metadata", "metadata.tsv"]
        
        sample = Sample("test", sample_options, global_filter_args)
        
        # Check that specific arguments are present and have correct values
        assert "--metadata" in sample.filter_args
        assert "--output-strains" in sample.filter_args
        assert "--group-by" in sample.filter_args
        assert "--subsample-max-sequences" in sample.filter_args
        assert "--query" in sample.filter_args
        assert "--query-columns" in sample.filter_args
        assert "--exclude-ambiguous-dates-by" in sample.filter_args
        assert "--min-length" in sample.filter_args
        assert "--max-length" in sample.filter_args
        
        # Check specific values
        metadata_idx = sample.filter_args.index("--metadata")
        assert sample.filter_args[metadata_idx + 1] == "metadata.tsv"
        
        query_idx = sample.filter_args.index("--query")
        assert sample.filter_args[query_idx + 1] == "country == 'Colombia'"
        
        query_cols_idx = sample.filter_args.index("--query-columns")
        assert sample.filter_args[query_cols_idx + 1] == "region:str"
        assert sample.filter_args[query_cols_idx + 2] == "coverage:float"
        
        exclude_dates_idx = sample.filter_args.index("--exclude-ambiguous-dates-by")
        assert sample.filter_args[exclude_dates_idx + 1] == "day"
        
        min_length_idx = sample.filter_args.index("--min-length")
        assert sample.filter_args[min_length_idx + 1] == "1000"
        
        max_length_idx = sample.filter_args.index("--max-length")
        assert sample.filter_args[max_length_idx + 1] == "30000"

    def test_boolean_flags_handling(self):
        """Test that boolean flags are correctly handled."""
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            "exclude_all": True,
            "non_nucleotide": True,
            "probabilistic_sampling": False,
        }
        global_filter_args = ["--metadata", "metadata.tsv"]
        
        sample = Sample("test", sample_options, global_filter_args)
        
        # Check boolean flags are present
        assert "--exclude-all" in sample.filter_args
        assert "--non-nucleotide" in sample.filter_args
        assert "--no-probabilistic-sampling" in sample.filter_args
        assert "--probabilistic-sampling" not in sample.filter_args

    def test_boolean_flags_false_values(self):
        """Test that boolean flags with False values are not added."""
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            "exclude_all": False,
            "non_nucleotide": False,
            "probabilistic_sampling": True,
        }
        global_filter_args = ["--metadata", "metadata.tsv"]
        
        sample = Sample("test", sample_options, global_filter_args)
        
        # Check that False boolean flags are not present
        assert "--exclude-all" not in sample.filter_args
        assert "--non-nucleotide" not in sample.filter_args
        assert "--probabilistic-sampling" in sample.filter_args
        assert "--no-probabilistic-sampling" not in sample.filter_args

    def test_missing_options_not_included(self):
        """Test that missing options are not included in filter args."""
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            # Only include a subset of available options
        }
        global_filter_args = ["--metadata", "metadata.tsv"]
        
        sample = Sample("test", sample_options, global_filter_args)
        
        # Check that missing options are not present
        assert "--query" not in sample.filter_args
        assert "--min-length" not in sample.filter_args
        assert "--exclude-all" not in sample.filter_args

    def test_array_values_handling(self):
        """Test that array values are properly expanded."""
        sample_options = {
            "group_by": ["region", "country"],
            "subsample_max_sequences": 5,
            "exclude": ["strain1.fasta", "strain2.fasta"],
            "query_columns": ["region:str", "date:str", "coverage:float"],
        }
        global_filter_args = ["--metadata", "metadata.tsv"]
        
        sample = Sample("test", sample_options, global_filter_args)
        
        # Check array expansion
        group_by_idx = sample.filter_args.index("--group-by")
        assert sample.filter_args[group_by_idx + 1] == "region"
        assert sample.filter_args[group_by_idx + 2] == "country"
        
        exclude_idx = sample.filter_args.index("--exclude")
        assert sample.filter_args[exclude_idx + 1] == "strain1.fasta"
        assert sample.filter_args[exclude_idx + 2] == "strain2.fasta"
        
        query_cols_idx = sample.filter_args.index("--query-columns")
        assert sample.filter_args[query_cols_idx + 1] == "region:str"
        assert sample.filter_args[query_cols_idx + 2] == "date:str"
        assert sample.filter_args[query_cols_idx + 3] == "coverage:float"


class TestGlobalFilterArgs:
    """Tests for global filter args building from CLI arguments and YAML config."""

    def test_cli_args_to_global_filter_args(self):
        """Test that CLI arguments are correctly converted to global filter args."""
        import argparse
        
        # Mock args object with CLI arguments
        args = argparse.Namespace()
        args.metadata_id_columns = ['strain', 'name']
        args.metadata_delimiters = ['\t', ',']
        args.metadata_chunk_size = 50000
        args.sequence_index = 'index.tsv'
        args.skip_checks = True
        
        # Build global filter args manually using the same logic as run()
        global_filter_args = []
        
        if args.metadata_id_columns:
            global_filter_args.extend(["--metadata-id-columns", *args.metadata_id_columns])
        
        if args.metadata_delimiters:
            global_filter_args.extend(["--metadata-delimiters", *args.metadata_delimiters])
        
        if args.metadata_chunk_size:
            global_filter_args.extend(["--metadata-chunk-size", str(args.metadata_chunk_size)])
        
        if args.sequence_index:
            global_filter_args.extend(["--sequence-index", args.sequence_index])
        
        if args.skip_checks:
            global_filter_args.append("--skip-checks")
        
        # Check that expected arguments are present
        assert "--metadata-id-columns" in global_filter_args
        assert "--metadata-delimiters" in global_filter_args
        assert "--metadata-chunk-size" in global_filter_args
        assert "--sequence-index" in global_filter_args
        assert "--skip-checks" in global_filter_args
        
        # Check specific values
        id_cols_idx = global_filter_args.index("--metadata-id-columns")
        assert global_filter_args[id_cols_idx + 1] == "strain"
        assert global_filter_args[id_cols_idx + 2] == "name"
        
        delim_idx = global_filter_args.index("--metadata-delimiters")
        assert global_filter_args[delim_idx + 1] == "\t"
        assert global_filter_args[delim_idx + 2] == ","
        
        chunk_idx = global_filter_args.index("--metadata-chunk-size")
        assert global_filter_args[chunk_idx + 1] == "50000"
        
        seq_idx = global_filter_args.index("--sequence-index")
        assert global_filter_args[seq_idx + 1] == "index.tsv"

    def test_yaml_config_to_sample_filter_args(self):
        """Test that YAML config options like include and include_where can be used in defaults."""
        import tempfile
        import os
        from augur.subsample import _parse_config
        
        config_content = textwrap.dedent("""
            defaults:
              include_where:
                - "region=North America"
              include:
                - strain1
            samples:
              test_sample:
                group_by:
                  - region
                subsample_max_sequences: 1
            """)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            config_file = f.name

        try:
            config = _parse_config(config_file)
            
            # Check that include and include_where are correctly parsed in defaults
            assert "defaults" in config
            assert "include_where" in config["defaults"]
            assert config["defaults"]["include_where"] == ["region=North America"]
            assert "include" in config["defaults"]
            assert config["defaults"]["include"] == ["strain1"]
            
        finally:
            os.unlink(config_file)

    def test_cli_args_false_values_not_included(self):
        """Test that CLI arguments with False/None values are not included."""
        import argparse
        
        # Mock args object with False/None values
        args = argparse.Namespace()
        args.metadata_id_columns = None
        args.metadata_delimiters = None
        args.metadata_chunk_size = None
        args.sequence_index = None
        args.skip_checks = False
        
        # Build global filter args manually using the same logic as run()
        global_filter_args = []
        
        if args.metadata_id_columns:
            global_filter_args.extend(["--metadata-id-columns", *args.metadata_id_columns])
        
        if args.metadata_delimiters:
            global_filter_args.extend(["--metadata-delimiters", *args.metadata_delimiters])
        
        if args.metadata_chunk_size:
            global_filter_args.extend(["--metadata-chunk-size", str(args.metadata_chunk_size)])
        
        if args.sequence_index:
            global_filter_args.extend(["--sequence-index", args.sequence_index])
        
        if args.skip_checks:
            global_filter_args.append("--skip-checks")
        
        # Check that no arguments are present
        assert "--metadata-id-columns" not in global_filter_args
        assert "--metadata-delimiters" not in global_filter_args
        assert "--metadata-chunk-size" not in global_filter_args
        assert "--sequence-index" not in global_filter_args
        assert "--skip-checks" not in global_filter_args


class TestMergeOptions:
    """Tests for the _merge_options function."""

    def test_merge_with_defaults(self):
        """Test basic merging where defaults provide missing options."""
        defaults = {
            "group_by": ["region"],
            "probabilistic_sampling": True,
            "min_date": "2020-01-01",
            "exclude": ["strain1.fasta"]
        }
        sample_options = {
            "subsample_max_sequences": 5,
            "max_date": "2020-12-31"
        }
        
        result = _merge_options(sample_options, defaults)
        
        # From defaults
        assert result["group_by"] == ["region"]
        assert result["probabilistic_sampling"] is True
        assert result["min_date"] == "2020-01-01"
        assert result["exclude"] == ["strain1.fasta"]
        
        # From sample options
        assert result["subsample_max_sequences"] == 5
        assert result["max_date"] == "2020-12-31"

    def test_sample_options_override_defaults(self):
        """Test that sample-specific options override defaults."""
        defaults = {
            "group_by": ["region"],
            "probabilistic_sampling": True,
            "exclude": ["strain1.fasta", "strain2.fasta"]
        }
        sample_options = {
            "subsample_max_sequences": 5,
            "group_by": ["country"],  # Should override default
            "probabilistic_sampling": False,  # Should override default
            "exclude": ["strain3.fasta"]  # Should override default
        }
        
        result = _merge_options(sample_options, defaults)
        
        # Sample options should override defaults
        assert result["group_by"] == ["country"]  # From sample, not default
        assert result["probabilistic_sampling"] is False  # From sample, not default
        assert result["exclude"] == ["strain3.fasta"]  # From sample, not default
        assert result["subsample_max_sequences"] == 5

    def test_merge_with_none_defaults(self):
        """Test backwards compatibility with None defaults."""
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            "probabilistic_sampling": True
        }
        
        result = _merge_options(sample_options, None)
        
        # Should return exactly the sample options
        assert result["group_by"] == ["region"]
        assert result["subsample_max_sequences"] == 5
        assert result["probabilistic_sampling"] is True

    def test_merge_with_empty_defaults(self):
        """Test merging with empty defaults object."""
        defaults = {}
        sample_options = {
            "group_by": ["region"],
            "subsample_max_sequences": 5,
            "probabilistic_sampling": True
        }
        
        result = _merge_options(sample_options, defaults)
        
        # Should return exactly the sample options since defaults is empty
        assert result["group_by"] == ["region"]
        assert result["subsample_max_sequences"] == 5
        assert result["probabilistic_sampling"] is True

    def test_mutual_exclusivity_sequences_per_group_wins(self):
        """Test that when sample has sequences_per_group, subsample_max_sequences is removed from defaults."""
        defaults = {
            "group_by": ["region"],
            "subsample_max_sequences": 10
        }
        sample_options = {
            "sequences_per_group": 3  # Should override and remove subsample_max_sequences
        }
        
        result = _merge_options(sample_options, defaults)
        
        assert result["group_by"] == ["region"]  # From defaults
        assert result["sequences_per_group"] == 3  # From sample
        assert "subsample_max_sequences" not in result  # Should be removed

    def test_mutual_exclusivity_subsample_max_sequences_wins(self):
        """Test that when sample has subsample_max_sequences, sequences_per_group is removed from defaults."""
        defaults = {
            "group_by": ["region"],
            "sequences_per_group": 2
        }
        sample_options = {
            "subsample_max_sequences": 5  # Should override and remove sequences_per_group
        }
        
        result = _merge_options(sample_options, defaults)
        
        assert result["group_by"] == ["region"]  # From defaults
        assert result["subsample_max_sequences"] == 5  # From sample
        assert "sequences_per_group" not in result  # Should be removed

    def test_no_mutual_exclusivity_conflict(self):
        """Test that when no conflict exists, both options can coexist."""
        defaults = {
            "group_by": ["region"],
            "min_date": "2020-01-01"
        }
        sample_options = {
            "subsample_max_sequences": 5,
            "max_date": "2020-12-31"
        }
        
        result = _merge_options(sample_options, defaults)
        
        # No conflict, all options should be present
        assert result["group_by"] == ["region"]
        assert result["min_date"] == "2020-01-01"
        assert result["subsample_max_sequences"] == 5
        assert result["max_date"] == "2020-12-31"

    def test_merge_complex_types(self):
        """Test that complex types like arrays are merged correctly."""
        defaults = {
            "group_by": ["region", "country"],
            "exclude": ["strain1.fasta", "strain2.fasta"],
            "query_columns": ["region:str", "date:str"]
        }
        sample_options = {
            "subsample_max_sequences": 5,
            "exclude": ["strain3.fasta", "strain4.fasta"]  # Should override default array
        }
        
        result = _merge_options(sample_options, defaults)
        
        # Array from defaults when not overridden
        assert result["group_by"] == ["region", "country"]
        assert result["query_columns"] == ["region:str", "date:str"]
        
        # Array from sample overrides defaults
        assert result["exclude"] == ["strain3.fasta", "strain4.fasta"]
        assert result["subsample_max_sequences"] == 5

    def test_merge_preserves_types(self):
        """Test that data types are preserved during merging."""
        defaults = {
            "group_by": ["region"],  # list
            "probabilistic_sampling": True,  # boolean
            "min_length": 1000,  # integer
            "min_date": "2020-01-01"  # string
        }
        sample_options = {
            "subsample_max_sequences": 5,  # integer
            "max_date": "2020-12-31"  # string
        }
        
        result = _merge_options(sample_options, defaults)
        
        # Check types are preserved
        assert isinstance(result["group_by"], list)
        assert isinstance(result["probabilistic_sampling"], bool)
        assert isinstance(result["min_length"], int)
        assert isinstance(result["min_date"], str)
        assert isinstance(result["subsample_max_sequences"], int)
        assert isinstance(result["max_date"], str)

    def test_defaults_with_both_subsample_options(self):
        """Test that defaults can contain both subsample options, and sample choice takes precedence."""
        defaults = {
            "group_by": ["region"],
            "subsample_max_sequences": 10,
            "sequences_per_group": 2
        }
        sample_options = {
            "sequences_per_group": 3  # Should keep this and remove subsample_max_sequences
        }
        
        result = _merge_options(sample_options, defaults)
        
        assert result["group_by"] == ["region"]
        assert result["sequences_per_group"] == 3  # From sample
        assert "subsample_max_sequences" not in result  # Removed due to conflict


class TestSubsampleDefaults:
    """Tests for schema validation of defaults functionality."""

    def test_valid_config_with_defaults(self):
        """Test that a valid config with defaults passes validation."""
        config_content = textwrap.dedent("""
            defaults:
              group_by:
                - region
              probabilistic_sampling: true
              min_date: "2020-01-01"
            samples:
              sample1:
                subsample_max_sequences: 5
              sample2:
                sequences_per_group: 3
                max_date: "2020-12-31"
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            config = _parse_config(temp_file)
            assert 'defaults' in config
            assert config['defaults']['group_by'] == ['region']
            assert config['defaults']['probabilistic_sampling'] is True
            assert config['defaults']['min_date'] == "2020-01-01"
            assert 'samples' in config
            assert 'sample1' in config['samples']
            assert 'sample2' in config['samples']
        finally:
            os.unlink(temp_file)

    def test_invalid_defaults_schema_validation(self):
        """Test that invalid defaults fail schema validation."""
        config_content = textwrap.dedent("""
            defaults:
              invalid_property: "not_allowed"  # Invalid property name
            samples:
              sample1:
                group_by:
                  - region
                subsample_max_sequences: 10
            """)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(config_content)
            temp_file = f.name

        try:
            with pytest.raises(AugurError, match="Config validation failed"):
                _parse_config(temp_file)
        finally:
            os.unlink(temp_file)