import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtRemoteServiceFactory;
/*
 * Copyright 1999,2004 The Apache Software Foundation.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

public class EntryRetrievalExample {

    public static void main(String[] args) {

        UniProtRemoteServiceFactory factory = new UniProtRemoteServiceFactory();

        //Create entry retrieval service
        EntryRetrievalService entryRetrievalService = factory.getEntryRetrievalService();

        //Retrieve UniProt entry by its accession number
        UniProtEntry entry = (UniProtEntry) entryRetrievalService.getUniProtEntry("P00634");

        System.out.println("entry = " + entry);

        //If entry with a given accession number is not found, entry will be equal null
        if (entry != null) {
            System.out.println("entry = " + entry.getUniProtId().getValue());
            System.out.println("entry.getDescription().getProteinNames() = " + entry.getDescription().getProteinNames());
            System.out.println("entry.getDescription().getECNumbers() = " + entry.getDescription().getECNumbers());

        }

    }
}
