python3 mr_link_2_standalone.py \
            --reference_bed example_files/reference_cohort \
            --sumstats_exposure example_files/non_causal_exposure.txt \
            --sumstats_outcome example_files/non_causal_outcome.txt \
            --out example_of_a_non_causal_effect.txt


python3 mr_link_2_standalone.py \
            --reference_bed example_files/reference_cohort \
            --sumstats_exposure example_files/yes_causal_exposure.txt \
            --sumstats_outcome example_files/yes_causal_outcome.txt \
            --out example_of_a_causal_effect.txt


python3 mr_link_2_standalone.py \
            --reference_bed example_files/reference_cohort \
            --sumstats_exposure example_files/yes_causal_exposure.txt \
            --sumstats_outcome example_files/yes_causal_outcome.txt \
            --out example_of_a_causal_effect_regions_prespecified.txt\
            --prespecified_regions 2:101532661-102532661,2:101532661-103480976

