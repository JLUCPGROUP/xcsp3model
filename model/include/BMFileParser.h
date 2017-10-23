/*
 * BMFileParser.h
 *
 *  Created on: 2016年9月12日
 *      Author: leezear
 */

#ifndef BMFILEPARSER_H_
#define BMFILEPARSER_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>

#include "libxml/SAX2.h"
#include "libxml/xmlstring.h"

#define _EXPORTING
#ifdef _EXPORTING
#define API_DECLSPEC __declspec(dllexport)
#else
#define API_DECLSPEC __declspec(dllimport)
#endif

using namespace std;


string _bm_path;
bool hasfile;

void startElement(void *user_data, const xmlChar *name, const xmlChar **attrs) {
	if (!xmlStrcmp(name, (const xmlChar *) "BMFile")) {
		hasfile = true;
	}
}

void endElement(void *parser, const xmlChar *name) {
	if (!xmlStrcmp(name, BAD_CAST "BMFile")) {
		hasfile = false;
	}
}

static void OnCharacters(void *ctx, const xmlChar *ch, int len) {
	if (hasfile) {
		_bm_path = string(ch, ch + len);
	}
}

bool API_DECLSPEC FindBMPath(const string file_name) {
	errno_t err;
	FILE *f;
	err = fopen_s(&f, file_name.c_str(), "rb");
	char chars[1024];
	int res = fread(chars, 1, 4, f);
	if (res <= 0) {
		return false;
	}

	xmlInitParser();

	xmlSAXHandler SAXHander;
	memset(&SAXHander, 0, sizeof(xmlSAXHandler));
	SAXHander.initialized = XML_SAX2_MAGIC;
	SAXHander.startElement = startElement;
	SAXHander.endElement = endElement;
	SAXHander.characters = OnCharacters;

	xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(&SAXHander, NULL, chars,
			res, NULL);

	while ((res = fread(chars, 1, sizeof(chars) - 1, f)) > 0) {
		if (xmlParseChunk(ctxt, chars, res, 0)) {
			xmlParserError(ctxt, "xmlParseChunk");
			return false;
		}
	}
	xmlParseChunk(ctxt, chars, 0, 1);

	xmlFreeParserCtxt(ctxt);
	xmlCleanupParser();
	fclose(f);
	return true;
}
#endif /* BMFILEPARSER_H_ */
