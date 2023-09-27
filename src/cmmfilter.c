/*
   gcc xmltest.c $(xml2-config --cflags) $(xml2-config --libs)

   Filter to select a specific chromosome from a cmm file. Possibly doing something else as well.
   It is probably a much better idea to regenerate a new cmm-file directly from coordinates and
   extra data.

*/

#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>

  int
main(int argc, char **argv)
{
  xmlDoc         *document;
  xmlNode        *root, *first_child, *node, *temp;
  char           *filename, *outfilename;

  if (argc < 3) {
    fprintf(stderr, "Usage: %s in.cmm out.cmm\n", argv[0]);
    return 1;
  }

  filename = argv[1];
  outfilename = argv[2];


  document = xmlReadFile(filename, NULL, 0);
  root = xmlDocGetRootElement(document);
  fprintf(stdout, "Root is <%s> (%i)\n", root->name, root->type);
  first_child = root->children;
  size_t nodeNumber = 0;

  int maxId = 2489;
  int minId = 0;

  size_t nFreeMarkers = 0;
  size_t nMarkers = 0;
  size_t nLinks = 0;

  for (node = first_child; node; node = node->next) {
    //        fprintf(stdout, "\t Child is <%s> (%i)\n", node->name, node->type);
    int remove = 0;

    if(node->type == 1) // Link
    {
      nMarkers++;
      // xmlChar * x = xmlGetProp(node, (const unsigned char *) "radius");
      //        printf("radius=%s\n", (char *) x);
      xmlChar * xid = xmlGetProp(node, (const unsigned char *) "id");
      xmlChar * xid1 = xmlGetProp(node, (const unsigned char *) "id1");
      xmlChar * xid2 = xmlGetProp(node, (const unsigned char *) "id2");

      if(xid != NULL) // A marker has an 'id'
      {
        int id = atoi((char*) xid);
        if(id < minId)
        {
          remove = 1;
        }
        if(id > maxId)
        {
          remove = 1;
        }
        if(remove)
        {
          nFreeMarkers++;
        }
      }

      if(xid1 != NULL) // A link has 'id1' and 'id2'
      {
        nLinks++;
        int id1 = atoi((char*) xid1);
        int id2 = atoi((char*) xid2);
        if(id1 < minId || id1 > maxId)
        {
          remove= 1;
        }
        if(id2 < minId || id2 > maxId)
        {
          remove = 1;
        }
      }

      nodeNumber++;
    }

    if(remove == 1)
    {
      temp = node;
      node = node->next;

      xmlUnlinkNode(temp);
      xmlFreeNode(temp);
    }
  }

  printf("Found: %zu markers, removed %zu.\n", nMarkers, nFreeMarkers);

  printf("Writing to %s\n", outfilename);
  xmlSaveFileEnc(outfilename, document, "UTF-8");
  fprintf(stdout, "...\n");
  return 0;
}
