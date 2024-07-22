#include "calc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> // isspace 함수를 사용하기 위해 추가

#define INITIAL_BUFFER_SIZE 1024

// 한 줄을 읽어 동적으로 할당된 메모리에 저장합니다.
char *read_line(FILE *file) {
    size_t buffer_size = INITIAL_BUFFER_SIZE;
    size_t length = 0;
    char *buffer = (char *)malloc(buffer_size);
    if (!buffer) {
        perror("메모리 할당 실패");
        return NULL;
    }

    while (fgets(buffer + length, buffer_size - length, file)) {
        length += strlen(buffer + length);
        // 줄 끝에 개행 문자가 있는지 확인합니다.
        if (buffer[length - 1] == '\n') {
            // 개행 문자를 널 문자로 대체하여 줄을 완성합니다.
            buffer[length - 1] = '\0';
            return buffer;
        }

        // 버퍼 크기를 두 배로 늘립니다.
        buffer_size *= 2;
        char *new_buffer = (char *)realloc(buffer, buffer_size);
        if (!new_buffer) {
            free(buffer);
            perror("메모리 재할당 실패");
            return NULL;
        }
        buffer = new_buffer;
    }

    // 파일의 끝에 도달했으나 줄이 비어있지 않은 경우 (마지막 줄 처리)
    if (length > 0) {
        return buffer;
    }

    // 줄이 비어있으면 메모리 해제
    free(buffer);
    return NULL;
}

// 파일을 읽어 각 줄을 동적 배열에 저장합니다.
char **read_file(const char *filename, int *num_lines) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("파일 열기 실패");
        return NULL;
    }

    // 초기 줄 수 용량 설정
    size_t capacity = 10; // 초기 용량을 적절히 설정할 수 있습니다.
    *num_lines = 0;
    char **lines = (char **)malloc(capacity * sizeof(char *));
    if (!lines) {
        perror("메모리 할당 실패");
        fclose(file);
        return NULL;
    }

    // 파일을 끝까지 읽어 각 줄을 배열에 저장
    char *line;
    while ((line = read_line(file)) != NULL) {
        // 용량이 부족할 경우 배열 크기를 늘립니다.
        if (*num_lines >= capacity) {
            capacity *= 2;
            char **new_lines = (char **)realloc(lines, capacity * sizeof(char *));
            if (!new_lines) {
                perror("메모리 재할당 실패");
                fclose(file);
                return NULL;
            }
            lines = new_lines;
        }

        // 줄의 마지막 공백 문자를 제거합니다.
        size_t len = strlen(line);
        while (len > 0 && isspace(line[len - 1])) {
            line[--len] = '\0';
        }

        lines[*num_lines] = line;
        (*num_lines)++;
    }

    fclose(file);
    return lines;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "사용법: %s <파일 경로>\n", argv[0]);
        return 1;
    }

    int num_lines;
    char **lines = read_file(argv[1], &num_lines);
    if (!lines) {
        return 1;
    }

    // 저장된 줄을 출력합니다.
    for (int i = 0; i < num_lines; i++) {
        printf("Line %d: %s\n", i + 1, lines[i]);
        //polys[i] = eval_string(lines[i]); use , 문자열 수정도 생각해볼것.
        //free(lines[i]); // 각 줄의 메모리 해제
    }
    
    Poly p;

    //char *ringfile = lines[0];
    init_calc(lines[0], 1);  //대신 create_ring( ) 쓸것 calc.c 651line 참조
    show_ring(CurrentRing);

    /*
    parse_string = lines[1];
    parse_string_index = 0;
    yyparse();
    print_poly(result);
    */
    //while ( 1 ) {
    //    fgets(buf,sizeof(buf),stdin);
      //  p = eval_string(lines[1]);
      //  print_poly(p); printf("\n");
    //}
    //print_poly(result);
    //free(lines); // 배열의 메모리 해제

    return 0;
}

    




/*
int main(int argc,char **argv)
{
  char *ringfile;
  int from_string;

  if ( argc == 1 )
    ringfile = 0;
  else if ( !strcmp(argv[1],"-s") ) {
    from_string = 1;
    ringfile = argv[2];
  } else {
    from_string = 0;
    ringfile = argv[1];
  }
  init_calc(ringfile,from_string);
  show_ring(CurrentRing);
  Input = stdin;
  while ( 1 ) {
    parse_string = 0;
    yyparse();
    print_poly(result); printf("\n");
  }
}
*/