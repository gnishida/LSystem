#include <glew.h>
#include "Shader.h"
#include <iostream>
#include <QFile>
#include <QTextStream>

uint Shader::createProgram(const char *vertex_file, const char* fragment_file)
{
	std::string source;
	loadTextFile(vertex_file, source);
	vertex_shader = loadShader(source, GL_VERTEX_SHADER);

	loadTextFile(fragment_file, source);
	fragment_shader = loadShader(source,GL_FRAGMENT_SHADER);

	// create program
	program = glCreateProgram();
	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glBindFragDataLocation(program, 0, "outputF");
	glLinkProgram(program);
	
	{
		int infologLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH,&infologLength);
		if (infologLength > 0){
			char *infoLog= (char *)malloc(infologLength);
			glGetProgramInfoLog(program, infologLength, NULL, infoLog);
			printf("%s\n",infoLog);
			free(infoLog);
		}
	}

	return program;
}

void Shader::loadTextFile(const QString& fileName, std::string& str)
{
	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly)) {
		printf("ERROR: loadTexFile: %s\n", file.errorString().toAscii().constData());
	}

	QTextStream in(&file);
	QString text;
	while(!in.atEnd()) {
		QString line = in.readLine();    
		text+=line+"\n"; 
	}
	file.close();
	str = std::string(text.toAscii().constData());
}

/**
 * Load a shader from the source.
 * Specify the corresponding mode, GL_VERTEX_SHADER or GL_FRAGMENT_SHADER.
 */
GLuint Shader::loadShader(std::string& source, GLuint mode)
{
	GLenum err;
	uint id = glCreateShader(mode);
	const char* csource = source.c_str();
	glShaderSource(id, 1, &csource, NULL);
	glCompileShader(id);

	int infologLength = 0;
	glGetShaderiv(id, GL_INFO_LOG_LENGTH, &infologLength);
	if (infologLength > 1) {
		char* infoLog = (char*)malloc(infologLength);
		glGetShaderInfoLog(id, 1000, NULL, infoLog);
		printf("Compile status %d:\n %s\n", infologLength, infoLog);
		free(infoLog);
		exit(1);
	} else {
		printf("Compile status: OK\n");
	}

	return id;
}

GLboolean Shader::printShaderInfoLog(GLuint shader, const char *str)
{
	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if (status == GL_FALSE) std::cerr << "Compile Error in " << str << std::endl;
  
	GLsizei bufSize;
	glGetShaderiv(shader, GL_INFO_LOG_LENGTH , &bufSize);
  
	if (bufSize > 1) {
		GLchar *infoLog = new GLchar[bufSize];
		GLsizei length;
		glGetShaderInfoLog(shader, bufSize, &length, infoLog);
		std::cerr << infoLog << std::endl;
		delete[] infoLog;
	}

	return (GLboolean)status;
}

GLboolean Shader::printProgramInfoLog(GLuint program)
{
	GLint status;
	glGetProgramiv(program, GL_LINK_STATUS, &status);
	if (status == GL_FALSE) std::cerr << "Link Error" << std::endl;
  
	GLsizei bufSize;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH , &bufSize);
	
	if (bufSize > 1) {
		GLchar *infoLog = new GLchar[bufSize];
		GLsizei length;
		glGetProgramInfoLog(program, bufSize, &length, infoLog);
		std::cerr << infoLog << std::endl;
		delete[] infoLog;
	}

	return (GLboolean)status;
}