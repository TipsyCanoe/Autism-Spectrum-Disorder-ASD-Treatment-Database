--
-- PostgreSQL database dump
--

-- Dumped from database version 17.4
-- Dumped by pg_dump version 17.2

-- Started on 2025-03-20 17:13:43

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET transaction_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- TOC entry 218 (class 1259 OID 24577)
-- Name: treatment_data; Type: TABLE; Schema: public; Owner: postgres
--

CREATE TABLE public.treatment_data (
    id integer NOT NULL,
    treatment character varying(255) NOT NULL,
    age character varying(50),
    gender character varying(50),
    duration character varying(255),
    n integer,
    primary_outcome_area character varying(255)
);


ALTER TABLE public.treatment_data OWNER TO postgres;

--
-- TOC entry 4899 (class 0 OID 0)
-- Dependencies: 218
-- Name: TABLE treatment_data; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON TABLE public.treatment_data IS 'Stores information about treatment methods and their applications';


--
-- TOC entry 4900 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.treatment; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.treatment IS 'Name of the treatment method';


--
-- TOC entry 4901 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.age; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.age IS 'Age range or mean age of participants';


--
-- TOC entry 4902 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.gender; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.gender IS 'Gender distribution of participants';


--
-- TOC entry 4903 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.duration; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.duration IS 'Duration of the treatment protocol';


--
-- TOC entry 4904 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.n; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.n IS 'Number of participants in the study';


--
-- TOC entry 4905 (class 0 OID 0)
-- Dependencies: 218
-- Name: COLUMN treatment_data.primary_outcome_area; Type: COMMENT; Schema: public; Owner: postgres
--

COMMENT ON COLUMN public.treatment_data.primary_outcome_area IS 'Primary area where outcomes were measured';


--
-- TOC entry 217 (class 1259 OID 24576)
-- Name: treatment_data_id_seq; Type: SEQUENCE; Schema: public; Owner: postgres
--

CREATE SEQUENCE public.treatment_data_id_seq
    AS integer
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER SEQUENCE public.treatment_data_id_seq OWNER TO postgres;

--
-- TOC entry 4906 (class 0 OID 0)
-- Dependencies: 217
-- Name: treatment_data_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: postgres
--

ALTER SEQUENCE public.treatment_data_id_seq OWNED BY public.treatment_data.id;


--
-- TOC entry 4742 (class 2604 OID 24580)
-- Name: treatment_data id; Type: DEFAULT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.treatment_data ALTER COLUMN id SET DEFAULT nextval('public.treatment_data_id_seq'::regclass);


--
-- TOC entry 4893 (class 0 OID 24577)
-- Dependencies: 218
-- Data for Name: treatment_data; Type: TABLE DATA; Schema: public; Owner: postgres
--

COPY public.treatment_data (id, treatment, age, gender, duration, n, primary_outcome_area) FROM stdin;
3	Arbaclofen (STX209)	5-21 y	124 M: 26 F	12 weeks	150	Social behaviors.
4	Arbaclofen (STX209)	6-39 y	55 M: 8 F	15 months	63	Irritability in fragile x syndrome (FXS)
1	Anodal transcranial direct current stimulation (tDCS)	5-8 y	100% M	8 weeks total (5 consecutive days tDCS, 1 week assessment, 4-week washout, 5 consecutive days days of tDCS)	20	N/A
2	Arbaclofen	5-21 y; mean: 11 y	124 M: 26 F	12 weeks	150	N/A
\.


--
-- TOC entry 4907 (class 0 OID 0)
-- Dependencies: 217
-- Name: treatment_data_id_seq; Type: SEQUENCE SET; Schema: public; Owner: postgres
--

SELECT pg_catalog.setval('public.treatment_data_id_seq', 4, true);


--
-- TOC entry 4746 (class 2606 OID 24584)
-- Name: treatment_data treatment_data_pkey; Type: CONSTRAINT; Schema: public; Owner: postgres
--

ALTER TABLE ONLY public.treatment_data
    ADD CONSTRAINT treatment_data_pkey PRIMARY KEY (id);


--
-- TOC entry 4743 (class 1259 OID 24586)
-- Name: idx_outcome_area; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_outcome_area ON public.treatment_data USING btree (primary_outcome_area);


--
-- TOC entry 4744 (class 1259 OID 24585)
-- Name: idx_treatment_name; Type: INDEX; Schema: public; Owner: postgres
--

CREATE INDEX idx_treatment_name ON public.treatment_data USING btree (treatment);


-- Completed on 2025-03-20 17:13:44

--
-- PostgreSQL database dump complete
--

