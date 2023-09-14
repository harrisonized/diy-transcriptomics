These are unaltered scripts from [https://diytranscriptomics.com](https://diytranscriptomics.com)

Some of the antipatterns you see in the scripts are listed here:

1. code duplication
2. hard-coding: paths, data, column names, etc.
3. missing references: some lines reference objects instantiated in previous scripts, requiring you to save all your RData when you close RStudio
4. lack of modularity: data and figures are all input and output into the same directory, making it difficult to figure out which files are used
5. manual intervention: instead of using `ggsave()` or `png() plot() dev.off()`, the script expects you to run each line in RStudio and use the [Export] button, making figure generation non-deterministic
6. improper object instantiation: filtered dataframes or data subsets are saved in permanent variables, signaling they might be important, but then they are never used again
7. lack of abstraction: no examples of writing reusable functions
8. improper scoping: importing entire tidyverse package prevents you from knowing which package your function comes from