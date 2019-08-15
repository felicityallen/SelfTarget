exporting data:
(based on https://medium.com/@pjausovec/exporting-mongodb-collections-to-json-808f3df6e988)

```bash
docker run -it -v $(PWD):/data mongo /bin/bash
# inside Docker
mongoexport --uri=mongodb://cellgeni:cellgeni@172.27.82.34:32556/wge -c w_g_e -o dump.json --jsonArray
```