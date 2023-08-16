function data = generate()
    % Extract correlation data from the dataset

    % Logic for Spearman correlation calculations from the original script
    % Spearman correlation comparing GP edit distance to phonological and orthographic
    GPprobToPhono(a, 1) = corr(phonemeTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToletter(a, 1) = corr(letterTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToSem(a, 1) = corr(semanticTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGraph(a, 1) = corr(graphemeTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToFreq_SUBTLEXUS(a, 1) = corr(Freq_SUBTLEXUS(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToOLD20sum(a, 1) = corr(OLD20sum(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGPfreq(a, 1) = corr(GPfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToPhonographic_N(a, 1) = corr(Phonographic_N(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToSem_N_D(a, 1) = corr(Sem_N_D(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToBigramF_Avg_C_Log(a, 1) = corr(BigramF_Avg_C_Log(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToLength(a, 1) = corr(Length(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToNSyll(a, 1) = corr(NSyll(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToBiphon_Prob(a, 1) = corr(Biphon_Prob(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToOR(a, 1) = corr(ORTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobTogp(a, 1) = corr(gpTri(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToORprob(a, 1) = corr(ORprob(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToGfreq(a, 1) = corr (Gfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');
    GPprobToPfreq(a, 1) = corr (Pfreq(tri_mask), GPprob(tri_mask), 'Type', 'Spearman');

    % Save or return the calculated data (This might need further refinement based on the library's structure)
    % data = ...;
end
