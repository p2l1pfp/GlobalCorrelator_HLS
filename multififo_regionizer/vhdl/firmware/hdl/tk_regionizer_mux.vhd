library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity tk_regionizer_mux is
    port(
            ap_clk : IN STD_LOGIC;
            ap_rst : IN STD_LOGIC;
            ap_start : IN STD_LOGIC;
            ap_done : OUT STD_LOGIC;
            ap_idle : OUT STD_LOGIC;
            ap_ready : OUT STD_LOGIC;
            newevent : IN STD_LOGIC;
            tracks_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_out       : OUT w64s(NTKSORTED-1 downto 0);
            newevent_out     : OUT STD_LOGIC
    );
end tk_regionizer_mux;

architecture Behavioral of tk_regionizer_mux is

    signal regionized:        w64s(NPFREGIONS-1 downto 0);
    signal regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal regionized_roll:   std_logic := '0';

    signal sorted_out :        particles(NPFREGIONS*NTKSORTED-1 downto 0);
    signal sorted_out_valid :  std_logic_vector(NPFREGIONS*NTKSORTED-1 downto 0) := (others => '0');
    signal sorted_out_roll :   std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');

    signal muxed_out :        particles(NTKSORTED-1 downto 0);
    signal muxed_out_valid :  std_logic_vector(NTKSORTED-1 downto 0) := (others => '0');
    signal muxed_out_roll :   std_logic := '0';

begin

    tk_regionizer : entity work.tk_regionizer 
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => ap_start,
                             newevent => newevent,
                             tracks_in_0_0_V => tracks_in_0_0_V,
                             tracks_in_0_1_V => tracks_in_0_1_V,
                             tracks_in_1_0_V => tracks_in_1_0_V,
                             tracks_in_1_1_V => tracks_in_1_1_V,
                             tracks_in_2_0_V => tracks_in_2_0_V,
                             tracks_in_2_1_V => tracks_in_2_1_V,
                             tracks_in_3_0_V => tracks_in_3_0_V,
                             tracks_in_3_1_V => tracks_in_3_1_V,
                             tracks_in_4_0_V => tracks_in_4_0_V,
                             tracks_in_4_1_V => tracks_in_4_1_V,
                             tracks_in_5_0_V => tracks_in_5_0_V,
                             tracks_in_5_1_V => tracks_in_5_1_V,
                             tracks_in_6_0_V => tracks_in_6_0_V,
                             tracks_in_6_1_V => tracks_in_6_1_V,
                             tracks_in_7_0_V => tracks_in_7_0_V,
                             tracks_in_7_1_V => tracks_in_7_1_V,
                             tracks_in_8_0_V => tracks_in_8_0_V,
                             tracks_in_8_1_V => tracks_in_8_1_V,
                             tracks_out_0_V => regionized(0),
                             tracks_out_1_V => regionized(1),
                             tracks_out_2_V => regionized(2),
                             tracks_out_3_V => regionized(3),
                             tracks_out_4_V => regionized(4),
                             tracks_out_5_V => regionized(5),
                             tracks_out_6_V => regionized(6),
                             tracks_out_7_V => regionized(7),
                             tracks_out_8_V => regionized(8),
                             tracks_out_valid_0 => regionized_valid(0),
                             tracks_out_valid_1 => regionized_valid(1),
                             tracks_out_valid_2 => regionized_valid(2),
                             tracks_out_valid_3 => regionized_valid(3),
                             tracks_out_valid_4 => regionized_valid(4),
                             tracks_out_valid_5 => regionized_valid(5),
                             tracks_out_valid_6 => regionized_valid(6),
                             tracks_out_valid_7 => regionized_valid(7),
                             tracks_out_valid_8 => regionized_valid(8),
                             newevent_out => regionized_roll);


    gen_sorters: for isort in NPFREGIONS-1 downto 0 generate
        reg_sorter : entity work.stream_sort
                            generic map(NITEMS => NTKSORTED)
                            port map(ap_clk => ap_clk,
                                d_in => w64_to_particle(regionized(isort)),
                                valid_in => regionized_valid(isort),
                                roll => regionized_roll,
                                d_out => sorted_out((isort+1)*NTKSORTED-1 downto isort*NTKSORTED),
                                valid_out => sorted_out_valid((isort+1)*NTKSORTED-1 downto isort*NTKSORTED),
                                roll_out => sorted_out_roll(isort)
                            );
        end generate gen_sorters;

    mux: entity work.region_mux
                            generic map(NREGIONS => NPFREGIONS, 
                                        NITEMS   => NTKSORTED,
                                        OUTII    => PFII)
                            port map(ap_clk => ap_clk,
                                roll => sorted_out_roll(0),
                                d_in => sorted_out,
                                valid_in => sorted_out_valid,
                                d_out => muxed_out,
                                valid_out => muxed_out_valid,
                                roll_out  => muxed_out_roll);

    format: process(ap_clk)
        begin
            if rising_edge(ap_clk) then
                for i in 0 to NTKSORTED-1 loop
                    if muxed_out_valid(i) = '1' then
                        tracks_out(i) <= particle_to_w64(muxed_out(i));
                    else
                        tracks_out(i) <= (others => '0');
                    end if;
                end loop;
                newevent_out <= muxed_out_roll;
            end if;
        end process format;


 end Behavioral;
