#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    FILE* file = fopen("src/user_inputs.in", "r");
    if (!file) {
        printf("Could not open file\n");
        return 1;
    }
    char buffer[1024];
    while(fgets(buffer, sizeof(buffer), file)) {
        char* hash = strchr(buffer, '#');
        if (hash) *hash = '\0';
        
        char* bookmark;
        char* vtoken = strtok_r(buffer, " \t\r\n", &bookmark);
        char* ktoken = strtok_r(NULL, " \t\r\n", &bookmark);
        
        if(!vtoken || !ktoken) continue;
        
        if(strstr(ktoken, "lk_correction_flag")) {
            printf("Found lk_correction_flag: value=%s\n", vtoken);
        } else {
            printf("Other token: k=%s v=%s\n", ktoken, vtoken);
        }
    }
    fclose(file);
    return 0;
}
